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
import numpy
import errno
import pandas
import os.path
import warnings
import triangle
from uuid import uuid4
from shutil import rmtree

class raster2TIN:
    """ 
    This class is useful for building the Badlands surface grid from a rectangular grid (DEM). This grid  
    is used to generate the irregular surface (TIN) om which the interactions between surface processes 
    and underlying sedimentary rocks will be computed.
    
    The purpose of the class is:
        1. to read and store the regular grid coordinates in numpy arrays format
        2. to define the edges of the grid computation by adding ghost vertices
        
    Parameters
    ----------
    string : inputfile
        This is a string containing the path to the regular grid file.
        
    string : outputDir
        This is a string containing the path to the output directory.
        
    integer : rank
        Rank of processor.
        
    string : delimiter
        The delimiter between columns from the regular grid file. The regular file contains
        coordinates of each nodes and is ordered by row from SW to NE corner. The file has no 
        header.
        Default: ' '
        
    variable: resRecFactor
        This integer gives the factor that will be used to define the resolution of the
        irregular grid edges 
                    >> TIN edges resolution = DEM resolution x resRecFactor
        Default: 1
        
    variable: areaDelFactor
        This integer gives the factor that will be used to define the averaged area of the
        irregular grid delaunay cells 
                    >> TIN cells resolution = areaDelFactor x (TIN edges resolution)^2
        Default: 1
    """
    
    def __init__(self, inputfile=None, outputDir=None, rank=0, delimiter=' ', resRecFactor=1, areaDelFactor=1):
        
        if inputfile==None:
            raise RuntimeError('DEM input file name must be defined to construct Badlands irregular grid.')
        if not os.path.isfile(inputfile):
            raise RuntimeError('The DEM input file name cannot be found in your path.')
        self.inputfile = inputfile
        
        if outputDir==None:
            raise RuntimeError('Output directory name must be defined to use Badlands code.')
        
        if rank == 0:
            temp_path = os.path.dirname(outputDir)+'/'+str(uuid4())
            try:
                os.renames(outputDir, temp_path)
            except OSError as exception:
                if exception.errno != errno.ENOENT:
                    raise
            else:
                rmtree(temp_path)
            os.mkdir(outputDir)
        
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
        minX = self.rectX.min()
        maxX = self.rectX.max()
        minY = self.rectY.min()
        maxY = self.rectY.max() 

        # Defines TIN variables
        self.resEdges = resDEM*self.resRecFactor
        self.resEdges = max(self.resEdges,resDEM)
        self.areaDel = (self.resEdges**2)*self.areaDelFactor
        self.areaDel = max(self.resEdges**2,self.areaDel)
        
        #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
        # Build edges
        #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
        # North South edges
        self.nx = int((maxX-minX)/self.resEdges+1)
        e_x = numpy.linspace(minX,maxX,self.nx)
        tmp1 = numpy.zeros(self.nx)
        tmp2 = numpy.zeros(self.nx)
        tmp1.fill(minY)
        tmp2.fill(maxY)
        south = numpy.column_stack((e_x,tmp1))
        north = numpy.column_stack((e_x,tmp2))
        
        # East West edges
        self.ny = int((maxY-minY)/self.resEdges+1)
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
        
        #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
        # Export DEM grid parameters as numpy arrays
        #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
        self.rnx = int((maxX-minX)/resDEM+1)
        self.rny = int((maxY-minY)/resDEM+1)
        self.regX = numpy.linspace(minX,maxX,self.rnx)
        self.regY = numpy.linspace(minY,maxY,self.rny)
        self.regZ = numpy.zeros((self.rnx,self.rny), dtype=float)
        p = 0
        for j in range(self.rny-1,-1,-1):
            for i in range(self.rnx):
                self.regZ[i,j] = self.rectZ[p]
                p += 1
        
        return
    
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
