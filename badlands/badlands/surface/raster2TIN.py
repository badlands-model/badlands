##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to the creation of triangular irregular network (TIN)
from raster type digital elevation model (DEM).

.. image:: img/tin_raster.png
   :scale: 90 %
   :alt: TIN to raster
   :align: center

"""
import os
import glob
import h5py
import numpy
import errno
import pandas
import os.path
import warnings
import tribad as triangle
from shutil import rmtree
from scipy import interpolate
from scipy.spatial import cKDTree


class raster2TIN:
    """
    Class to build **badlands** surface grid from a rectangular grid (DEM).

    Note:
        This grid is used to generate the irregular surface (TIN) on which the interactions
        between surface processes and underlying sedimentary rocks will be computed.

    The purpose of the class is twofold:

    1. reading and storing the regular grid coordinates in numpy arrays format
    2. defining the **edges** of the grid computation by adding **ghost** nodes

    The regular file contains coordinates of each nodes and is ordered by row from **SW** to **NE**
    corner.

    Args:
        inputfile : (str) this is a string containing the path to the regular grid file.
        delimiter : (str) delimiter between columns from the regular grid file (default: ' ')
        resRecFactor : this integer gives the factor that will be used to define the resolution of the irregular grid edges (default: 1)
        areaDelFactor : This integer gives the factor that will be used to define the averaged area of the irregular grid delaunay cells (default: 1)

    Caution:
        * The input DEM file has no header.
        * The TIN edges resolution is equal to the initial DEM resolution times the chosen resRecFactor from the XML file.
        * The TIN cells resolution is equal to the areaDelFactor times the square of the TIN edges resolution.

    """

    def __init__(
        self, inputfile=None, delimiter=r"\s+", resRecFactor=1, areaDelFactor=1
    ):

        if inputfile == None:
            raise RuntimeError(
                "DEM input file name must be defined to construct Badlands irregular grid."
            )
        if not os.path.isfile(inputfile):
            raise RuntimeError("The DEM input file name cannot be found in your path.")
        self.inputfile = inputfile

        self.delimiter = delimiter

        if resRecFactor < 1:
            raise ValueError("TIN edges resolution factor needs to be at least 1.")
        self.resRecFactor = resRecFactor

        if areaDelFactor < 1:
            raise ValueError("TIN cell area factor needs to be at least 1.")
        self.areaDelFactor = areaDelFactor

        # Define class parameters
        self.nx = None
        self.ny = None
        self.areaDel = None
        self.resEdges = None
        self.rectX = None
        self.rectY = None
        self.rectZ = None
        self.edges = None
        self.edgesPt = None
        self.bounds = None
        self.boundsPt = None
        self.rnx = None
        self.rny = None
        self.regX = None
        self.regY = None
        self.regZ = None
        self.tinMesh = None
        self.bmask = None
        self.resdx = None

        # TIN creation
        self._triangulate_raster_from_file()

    def _raster_edges(self):
        """
        Using Pandas library to read the DEM file and allocating nodes and edges.
        This function also sets the TIN parameters.
        """

        # Read DEM file
        data = pandas.read_csv(
            self.inputfile,
            sep=self.delimiter,
            engine="c",
            header=None,
            na_filter=False,
            dtype=numpy.float,
            low_memory=False,
        )
        self.rectX = data.values[:, 0]
        self.rectY = data.values[:, 1]
        self.rectZ = data.values[:, 2]
        resDEM = self.rectX[1] - self.rectX[0]
        self.resdx = resDEM
        minX = self.rectX.min()
        maxX = self.rectX.max()
        minY = self.rectY.min()
        maxY = self.rectY.max()

        # Defines TIN variables
        self.resEdges = resDEM * self.resRecFactor
        self.resEdges = max(self.resEdges, resDEM)
        self.areaDel = (self.resEdges ** 2) * self.areaDelFactor
        self.areaDel = max(self.resEdges ** 2, self.areaDel)

        # North South edges
        self.nx = int(round((maxX - minX) / self.resEdges + 1))
        e_x = numpy.linspace(minX, maxX, self.nx)
        tmp1 = numpy.zeros(self.nx)
        tmp2 = numpy.zeros(self.nx)
        tmp1.fill(minY)
        tmp2.fill(maxY)
        south = numpy.column_stack((e_x, tmp1))
        north = numpy.column_stack((e_x, tmp2))

        # East West edges
        self.ny = int(round((maxY - minY) / self.resEdges + 1))
        e_y = numpy.linspace(minY + self.resEdges, maxY - self.resEdges, self.ny - 2)
        tmp1 = numpy.zeros(self.ny - 2)
        tmp2 = numpy.zeros(self.ny - 2)
        tmp1.fill(minX)
        tmp2.fill(maxX)
        east = numpy.column_stack((tmp1, e_y))
        west = numpy.column_stack((tmp2, e_y))

        # Merge edges together
        edges = []
        edges = numpy.vstack((south, east))
        edges = numpy.vstack((edges, north))
        edges = numpy.vstack((edges, west))
        self.edges = edges

        self.edgesPt = len(self.edges)

        self.rnx = int(round((maxX - minX) / resDEM + 1))
        self.rny = int(round((maxY - minY) / resDEM + 1))
        self.regX = numpy.linspace(minX, maxX, self.rnx)
        self.regY = numpy.linspace(minY, maxY, self.rny)
        self.regZ = numpy.reshape(self.rectZ, (self.rnx, self.rny), order="F")

    def _TIN_ghosts_bounds(self):
        """
        Extend the boundary of the TIN grid using ghost cells on the edges of the
        regular grid.
        """

        minX = self.rectX.min()
        maxX = self.rectX.max()
        minY = self.rectY.min()
        maxY = self.rectY.max()

        # North South bounds
        b_x = numpy.linspace(minX - self.resEdges, maxX + self.resEdges, self.nx + 2)
        tmp1 = numpy.zeros(self.nx + 2)
        tmp2 = numpy.zeros(self.nx + 2)
        tmp1.fill(minY - self.resEdges)
        tmp2.fill(maxY + self.resEdges)
        south = numpy.column_stack((b_x, tmp1))
        north = numpy.column_stack((b_x, tmp2))

        # East West edges
        b_y = numpy.linspace(minY, maxY, self.ny)
        tmp1 = numpy.zeros(self.ny)
        tmp2 = numpy.zeros(self.ny)
        tmp1.fill(minX - self.resEdges)
        tmp2.fill(maxX + self.resEdges)
        east = numpy.column_stack((tmp1, b_y))
        west = numpy.column_stack((tmp2, b_y))

        # Merge edges together
        bounds = []
        bounds = numpy.vstack((south, east))
        bounds = numpy.vstack((bounds, north))
        bounds = numpy.vstack((bounds, west))
        self.bounds = bounds

        self.boundsPt = len(self.bounds)

        return

    def _triangulate_raster_from_file(self):
        """
        Main function used to create the irregular TIN surface for Badlands based on
        a regular DEM file.
        """

        # Compute DEM edges
        self._raster_edges()

        # Compute TIN grid boundaries
        self._TIN_ghosts_bounds()

        # Create TIN
        tinPts = numpy.vstack((self.bounds, self.edges))
        self.tinMesh = triangle.triangulate(
            {"vertices": tinPts}, "Dqa" + str(self.areaDel)
        )
        ptsTIN = self.tinMesh["vertices"]

        # Check extent
        checkXmin = ptsTIN[self.boundsPt :, 0].min()
        checkXmax = ptsTIN[self.boundsPt :, 0].max()
        checkYmin = ptsTIN[self.boundsPt :, 1].min()
        checkYmax = ptsTIN[self.boundsPt :, 1].max()

        if checkXmin < self.rectX.min() or checkXmax > self.rectX.max():
            raise ValueError(
                "Error in defining the X boundary nodes, you will need to adjust the resolution value"
            )
        if checkYmin < self.rectY.min() or checkYmax > self.rectY.max():
            raise ValueError(
                "Error in defining the Y boundary nodes, you will need to adjust the resolution value"
            )

        # Add masks
        self.bmask = numpy.zeros(len(ptsTIN[:, 0]))
        self.bmask[: self.boundsPt] = 1

        return

    def load_hdf5(self, restartFolder, timestep, tXY):
        """
        This function reads and allocates **badlands** parameters & variables from a previous
        simulation. These parameters are obtained from a specific time HDF5 file.

        Important:
            This function is only called when a **restart simulation** is ran... |:bomb:|

        Args:
            restartFolder: restart folder name as defined in tge XML input file.
            timestep: time step to load defined based on output interval number.
            tXY: 2D numpy array TIN grid local coordinates.

        Returns
        -------
        elev
            numpy 1D array containing the updated elevation from the restart model.
        cum
            numpy 1D array containing the updated erosion/deposition values from the restart model.
        """

        if os.path.exists(restartFolder):
            folder = restartFolder + "/h5/"
            fileCPU = "tin.time%s.hdf5" % timestep
            restartncpus = len(glob.glob1(folder, fileCPU))
            if restartncpus == 0:
                raise ValueError(
                    "The requested time step for the restart simulation cannot be found in the restart folder."
                )
        else:
            raise ValueError(
                "The restart folder is missing or the given path is incorrect."
            )

        df = h5py.File("%s/h5/tin.time%s.hdf5" % (restartFolder, timestep), "r")
        coords = numpy.array((df["/coords"]))
        cumdiff = numpy.array((df["/cumdiff"]))
        cumhill = numpy.array((df["/cumhill"]))
        cumfail = numpy.array((df["/cumfail"]))
        x, y, z = numpy.hsplit(coords, 3)
        c = cumdiff
        h = cumhill
        f = cumfail

        XY = numpy.column_stack((x, y))
        tree = cKDTree(XY)
        distances, indices = tree.query(tXY, k=3)

        if len(z[indices].shape) == 3:
            z_vals = z[indices][:, :, 0]
            c_vals = c[indices][:, :, 0]
            h_vals = h[indices][:, :, 0]
            f_vals = f[indices][:, :, 0]
        else:
            z_vals = z[indices]
            c_vals = c[indices]
            h_vals = h[indices]
            f_vals = f[indices]

        with numpy.errstate(divide="ignore"):
            elev = numpy.average(z_vals, weights=(1.0 / distances), axis=1)
            cum = numpy.average(c_vals, weights=(1.0 / distances), axis=1)
            hcum = numpy.average(h_vals, weights=(1.0 / distances), axis=1)
            fcum = numpy.average(f_vals, weights=(1.0 / distances), axis=1)

        onIDs = numpy.where(distances[:, 0] == 0)[0]
        if len(onIDs) > 0:
            if len(z[indices].shape) == 3:
                elev[onIDs] = z[indices[onIDs, 0], 0]
                cum[onIDs] = c[indices[onIDs, 0], 0]
                hcum[onIDs] = h[indices[onIDs, 0], 0]
                fcum[onIDs] = f[indices[onIDs, 0], 0]
            else:
                elev[onIDs] = z[indices[onIDs, 0]]
                cum[onIDs] = c[indices[onIDs, 0]]
                hcum[onIDs] = h[indices[onIDs, 0]]
                fcum[onIDs] = f[indices[onIDs, 0]]

        return elev, cum, hcum, fcum

    def load_hdf5_flex(self, restartFolder, timestep, tXY):
        """
        This function reads and allocates **badlands** parameters & variables from a previous
        simulation with flexural isostasy turned on.

        These parameters are obtained from a specific time HDF5 file.

        Args:
            restartFolder: restart folder name as defined in tge XML input file.
            timestep: time step to load defined based on output interval number.
            tXY: 2D numpy array TIN grid local coordinates.

        Returns
        -------
        elev
            numpy 1D array containing the updated elevation from the restart model.
        cum
            numpy 1D array containing the updated erosion/deposition values from the restart model.
        fcum
            numpy 1D array containing the cumulative flexural isostasy values from the restart model.


        Note:
            Flexural isostasy is obtained from the gFlex_ package available on Github!

            Wickert, A. D. (2016), Open-source modular solutions for flexural isostasy: gFlex v1.0,
            Geosci. Model Dev., 9(3), 997–1017, `doi:10.5194/gmd-9-997-2016`_.

        .. _`doi:10.5194/gmd-9-997-2016`:  https://doi.org/10.5194/gmd-9-997-2016
        .. _gflex: https://github.com/awickert/gFlex

        """

        if os.path.exists(restartFolder):
            folder = restartFolder + "/h5/"
            fileCPU = "tin.time%s.hdf5" % timestep
            restartncpus = len(glob.glob1(folder, fileCPU))
            if restartncpus == 0:
                raise ValueError(
                    "The requested time step for the restart simulation cannot be found in the restart folder."
                )
        else:
            raise ValueError(
                "The restart folder is missing or the given path is incorrect."
            )

        df = h5py.File("%s/h5/tin.time%s.hdf5" % (restartFolder, timestep), "r")
        coords = numpy.array((df["/coords"]))
        cumdiff = numpy.array((df["/cumdiff"]))
        cumhill = numpy.array((df["/cumhill"]))
        cumfail = numpy.array((df["/cumfail"]))
        cumflex = numpy.array((df["/cumflex"]))
        x, y, z = numpy.hsplit(coords, 3)
        c = cumdiff
        h = cumhill
        s = cumfail
        f = cumflex

        XY = numpy.column_stack((x, y))
        tree = cKDTree(XY)
        distances, indices = tree.query(tXY, k=3)

        if len(z[indices].shape) == 3:
            z_vals = z[indices][:, :, 0]
            c_vals = c[indices][:, :, 0]
            h_vals = h[indices][:, :, 0]
            s_vals = s[indices][:, :, 0]
            f_vals = f[indices][:, :, 0]
        else:
            z_vals = z[indices]
            c_vals = c[indices]
            s_vals = s[indices]
            h_vals = h[indices]
            f_vals = f[indices]

        distances[distances < 0.0001] = 0.0001
        with numpy.errstate(divide="ignore"):
            elev = numpy.average(z_vals, weights=(1.0 / distances), axis=1)
            cum = numpy.average(c_vals, weights=(1.0 / distances), axis=1)
            hcum = numpy.average(h_vals, weights=(1.0 / distances), axis=1)
            scum = numpy.average(s_vals, weights=(1.0 / distances), axis=1)
            fcum = numpy.average(f_vals, weights=(1.0 / distances), axis=1)

        onIDs = numpy.where(distances[:, 0] <= 0.0001)[0]
        if len(onIDs) > 0:
            if len(z[indices].shape) == 3:
                elev[onIDs] = z[indices[onIDs, 0], 0]
                cum[onIDs] = c[indices[onIDs, 0], 0]
                hcum[onIDs] = h[indices[onIDs, 0], 0]
                scum[onIDs] = s[indices[onIDs, 0], 0]
                fcum[onIDs] = f[indices[onIDs, 0], 0]
            else:
                elev[onIDs] = z[indices[onIDs, 0]]
                cum[onIDs] = c[indices[onIDs, 0]]
                hcum[onIDs] = h[indices[onIDs, 0]]
                scum[onIDs] = s[indices[onIDs, 0]]
                fcum[onIDs] = f[indices[onIDs, 0]]

        return elev, cum, hcum, scum, fcum
