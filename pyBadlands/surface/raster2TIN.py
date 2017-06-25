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
"""
import os
import glob
import h5py
import numpy
import errno
import pandas
import os.path
import warnings
import triangle
from uuid import uuid4
from shutil import rmtree
from scipy import interpolate
from scipy.spatial import cKDTree

class raster2TIN:
    """
    This class is useful for building the Badlands surface grid from a rectangular grid (DEM). This grid
    is used to generate the irregular surface (TIN) on which the interactions between surface processes
    and underlying sedimentary rocks will be computed.

    The purpose of the class is:
        1. to read and store the regular grid coordinates in numpy arrays format
        2. to define the edges of the grid computation by adding ghost vertices

    Parameters
    ----------
    inputfile : string
        This is a string containing the path to the regular grid file.

    rank : integer
        Rank of processor.

    delimiter : string
        The delimiter between columns from the regular grid file. The regular file contains
        coordinates of each nodes and is ordered by row from SW to NE corner. The file has no
        header.

        Default: ' '

    resRecFactor
        This integer gives the factor that will be used to define the resolution of the
        irregular grid edges

            >> TIN edges resolution = DEM resolution x resRecFactor

        Default: 1

    areaDelFactor
        This integer gives the factor that will be used to define the averaged area of the
        irregular grid delaunay cells

            >> TIN cells resolution = areaDelFactor x (TIN edges resolution)^2

        Default: 1
    """

    def __init__(self, inputfile=None, rank=0, delimiter=r'\s+', resRecFactor=1, areaDelFactor=1):
        if inputfile==None:
            raise RuntimeError('DEM input file name must be defined to construct Badlands irregular grid.')
        if not os.path.isfile(inputfile):
            raise RuntimeError('The DEM input file name cannot be found in your path.')
        self.inputfile = inputfile

        self.delimiter = delimiter

        if resRecFactor < 1:
            raise ValueError( "TIN edges resolution factor needs to be at least 1." )
        self.resRecFactor = resRecFactor

        if areaDelFactor < 1:
            raise ValueError( "TIN cell area factor needs to be at least 1." )
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
        data = pandas.read_csv(self.inputfile, sep=self.delimiter, engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)
        self.rectX = data.values[:,0]
        self.rectY = data.values[:,1]
        self.rectZ = data.values[:,2]
        resDEM = self.rectX[1]-self.rectX[0]
        self.resdx = resDEM
        minX = self.rectX.min()
        maxX = self.rectX.max()
        minY = self.rectY.min()
        maxY = self.rectY.max()

        # Defines TIN variables
        self.resEdges = resDEM*self.resRecFactor
        self.resEdges = max(self.resEdges,resDEM)
        self.areaDel = (self.resEdges**2)*self.areaDelFactor
        self.areaDel = max(self.resEdges**2,self.areaDel)

        # North South edges
        self.nx = int(round((maxX-minX)/self.resEdges+1))
        e_x = numpy.linspace(minX,maxX,self.nx)
        tmp1 = numpy.zeros(self.nx)
        tmp2 = numpy.zeros(self.nx)
        tmp1.fill(minY)
        tmp2.fill(maxY)
        south = numpy.column_stack((e_x,tmp1))
        north = numpy.column_stack((e_x,tmp2))

        # East West edges
        self.ny = int(round((maxY-minY)/self.resEdges+1))
        e_y = numpy.linspace(minY+self.resEdges,maxY-self.resEdges,self.ny-2)
        tmp1 = numpy.zeros(self.ny-2)
        tmp2 = numpy.zeros(self.ny-2)
        tmp1.fill(minX)
        tmp2.fill(maxX)
        east = numpy.column_stack((tmp1,e_y))
        west = numpy.column_stack((tmp2,e_y))

        # Merge edges together
        edges = []
        edges = numpy.vstack((south,east))
        edges = numpy.vstack((edges,north))
        edges = numpy.vstack((edges,west))
        self.edges = edges

        self.edgesPt = len(self.edges)

        self.rnx = int(round((maxX-minX)/resDEM+1))
        self.rny = int(round((maxY-minY)/resDEM+1))
        self.regX = numpy.linspace(minX,maxX,self.rnx)
        self.regY = numpy.linspace(minY,maxY,self.rny)
        self.regZ = numpy.reshape(self.rectZ,(self.rnx,self.rny),order='F')

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
        b_x = numpy.linspace(minX-self.resEdges,maxX+self.resEdges,self.nx+2)
        tmp1 = numpy.zeros(self.nx+2)
        tmp2 = numpy.zeros(self.nx+2)
        tmp1.fill(minY-self.resEdges)
        tmp2.fill(maxY+self.resEdges)
        south = numpy.column_stack((b_x,tmp1))
        north = numpy.column_stack((b_x,tmp2))

        # East West edges
        b_y = numpy.linspace(minY,maxY,self.ny)
        tmp1 = numpy.zeros(self.ny)
        tmp2 = numpy.zeros(self.ny)
        tmp1.fill(minX-self.resEdges)
        tmp2.fill(maxX+self.resEdges)
        east = numpy.column_stack((tmp1,b_y))
        west = numpy.column_stack((tmp2,b_y))

        # Merge edges together
        bounds = []
        bounds = numpy.vstack((south,east))
        bounds = numpy.vstack((bounds,north))
        bounds = numpy.vstack((bounds,west))
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
        tinPts = numpy.vstack(( self.bounds, self.edges))
        self.tinMesh = triangle.triangulate( dict(vertices=tinPts),'Dqa'+str(self.areaDel))
        ptsTIN = self.tinMesh['vertices']

        # Check extent
        checkXmin = ptsTIN[self.boundsPt:,0].min()
        checkXmax = ptsTIN[self.boundsPt:,0].max()
        checkYmin = ptsTIN[self.boundsPt:,1].min()
        checkYmax = ptsTIN[self.boundsPt:,1].max()

        if checkXmin < self.rectX.min() or checkXmax > self.rectX.max():
            raise ValueError('Error in defining the X boundary nodes, you will need to adjust the resolution value')
        if checkYmin < self.rectY.min() or checkYmax > self.rectY.max():
            raise ValueError('Error in defining the Y boundary nodes, you will need to adjust the resolution value')

        # Add masks
        self.bmask = numpy.zeros(len(ptsTIN[:,0]))
        self.bmask[:self.boundsPt] = 1

        return

    def load_hdf5(self, restartFolder, timestep, tXY):
        """
        Read the HDF5 file for a given time step.

        Parameters
        ----------
        restartFolder
            Restart folder name.

        timestep
            Time step to load.

        tXY
            TIN grid local coordinates.

        Returns
        -------
        elev
            Numpy array containing the updated elevation from the restart model.

        cum
            Numpy array containing the updated erosion/deposition values from the restart model.
        """

        if os.path.exists(restartFolder):
            folder = restartFolder+'/h5/'
            fileCPU = 'tin.time%s.p*.hdf5'%timestep
            restartncpus = len(glob.glob1(folder,fileCPU))
            if restartncpus == 0:
                raise ValueError('The requested time step for the restart simulation cannot be found in the restart folder.')
        else:
            raise ValueError('The restart folder is missing or the given path is incorrect.')


        for i in range(0, restartncpus):
            df = h5py.File('%s/h5/tin.time%s.p%s.hdf5'%(restartFolder, timestep, i), 'r')
            coords = numpy.array((df['/coords']))
            cumdiff = numpy.array((df['/cumdiff']))
            cumhill = numpy.array((df['/cumhill']))
            if i == 0:
                x, y, z = numpy.hsplit(coords, 3)
                c = cumdiff
                h = cumhill
            else:
                c = numpy.append(c, cumdiff)
                h = numpy.append(h, cumhill)
                x = numpy.append(x, coords[:,0])
                y = numpy.append(y, coords[:,1])
                z = numpy.append(z, coords[:,2])

        XY = numpy.column_stack((x,y))
        tree = cKDTree(XY)
        distances, indices = tree.query(tXY, k=3)

        if len(z[indices].shape) == 3:
            z_vals = z[indices][:,:,0]
            c_vals = c[indices][:,:,0]
            h_vals = h[indices][:,:,0]
        else:
            z_vals = z[indices]
            c_vals = c[indices]
            h_vals = h[indices]

        with numpy.errstate(divide='ignore'):
            elev = numpy.average(z_vals,weights=(1./distances), axis=1)
            cum = numpy.average(c_vals,weights=(1./distances), axis=1)
            hcum = numpy.average(h_vals,weights=(1./distances), axis=1)

        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            elev[onIDs] = z[indices[onIDs,0]]
            cum[onIDs] = c[indices[onIDs,0]]
            hcum[onIDs] = h[indices[onIDs,0]]

        return elev, cum, hcum

    def load_hdf5_flex(self, restartFolder, timestep, tXY):
        """
        Read the HDF5 file for a given time step when flexural isostasy is on.

        Parameters
        ----------
        restartFolder
            Restart folder name.

        timestep
            Time step to load.

        tXY
            TIN grid local coordinates.

        Returns
        -------
        elev
            Numpy array containing the updated elevation from the restart model.

        cum
            Numpy array containing the updated erosion/deposition values from the restart model.

        cumf
            Numpy array containing the cumulative flexural isostasy values from the restart model.
        """

        if os.path.exists(restartFolder):
            folder = restartFolder+'/h5/'
            fileCPU = 'tin.time%s.p*.hdf5'%timestep
            restartncpus = len(glob.glob1(folder,fileCPU))
            if restartncpus == 0:
                raise ValueError('The requested time step for the restart simulation cannot be found in the restart folder.')
        else:
            raise ValueError('The restart folder is missing or the given path is incorrect.')


        for i in range(0, restartncpus):
            df = h5py.File('%s/h5/tin.time%s.p%s.hdf5'%(restartFolder, timestep, i), 'r')
            coords = numpy.array((df['/coords']))
            cumdiff = numpy.array((df['/cumdiff']))
            cumhill = numpy.array((df['/cumhill']))
            cumflex = numpy.array((df['/cumflex']))
            if i == 0:
                x, y, z = numpy.hsplit(coords, 3)
                c = cumdiff
                h = cumhill
                f = cumflex
            else:
                c = numpy.append(c, cumdiff)
                h = numpy.append(h, cumhill)
                f = numpy.append(f, cumflex)
                x = numpy.append(x, coords[:,0])
                y = numpy.append(y, coords[:,1])
                z = numpy.append(z, coords[:,2])

        XY = numpy.column_stack((x,y))
        tree = cKDTree(XY)
        distances, indices = tree.query(tXY, k=3)

        if len(z[indices].shape) == 3:
            z_vals = z[indices][:,:,0]
            c_vals = c[indices][:,:,0]
            h_vals = h[indices][:,:,0]
            f_vals = f[indices][:,:,0]
        else:
            z_vals = z[indices]
            c_vals = c[indices]
            h_vals = h[indices]
            f_vals = f[indices]

        distances[distances<0.0001] = 0.0001
        with numpy.errstate(divide='ignore'):
            elev = numpy.average(z_vals,weights=(1./distances), axis=1)
            cum = numpy.average(c_vals,weights=(1./distances), axis=1)
            hcum = numpy.average(h_vals,weights=(1./distances), axis=1)
            cumf = numpy.average(f_vals,weights=(1./distances), axis=1)

        onIDs = numpy.where(distances[:,0] <= 0.0001)[0]
        if len(onIDs) > 0:
            elev[onIDs] = z[indices[onIDs,0]]
            cum[onIDs] = c[indices[onIDs,0]]
            cumf[onIDs] = f[indices[onIDs,0]]
            hcum[onIDs] = h[indices[onIDs,0]]

        return elev, cum, hcum, cumf
