##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
We use a finite volume (**FV**) approach to solve the continuity equations driving
sediment erosion, transport and deposition in **badlands**.

The set of functions below are defining the **FV** discretisation strategy.
"""

import time
import numpy
import warnings
import tribad as triangle

import os
if 'READTHEDOCS' not in os.environ:
    from badlands import fvframe

class FVmethod:
    """
    This class builds paramters required for the Finite Volume mesh algorithm.

    It builds for the TIN's the **dual Delaunay-Voronoi framework** and computes some of
    the geometrical characteristics of the numerical mesh such as the *voronoi cell area*,
    an ordered list of *voronoi nodes* or the *voronoi edges* lengths to cite a few...

    .. image:: img/voro.png
       :scale: 90 %
       :alt: TIN grid
       :align: center

    Args:
        nodes : numpy floating array of 2D coordinates of TIN's nodes position.
        cells : numpy integer array of 3 indices defining TIN's cells connectivity.
        edges : numpy integer array of 2 indices defining TIN's edges connectivity.
    """
    def __init__(self, nodes, cells, edges):

        self.node_coords = nodes
        self.edges = edges
        self.cells = cells
        self.control_volumes = None
        self.neighbours = None
        self.vor_edges = None
        self.edge_length = None
        self.fillH = None
        self.partIDs = None
        self.maxNgbh = None
        self.localIDs = None
        self.outPts = None
        self.outCells = None

    def _FV_utils(self, lGIDs, verbose=False):
        """
        This function constructs the Finite Volume discretisation for each local triangularised grid.

        Args:
            lGIDs: numpy integer-type array filled with the global vertex IDs for each local grid located within the partition (including those on the edges).
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).
        """

        # Build the voronoi diagram (dual of the delaunay)
        walltime = time.clock()
        Vor_pts, Vor_edges = triangle.voronoi(self.node_coords)
        if verbose:
            print(" - build the voronoi diagram ", time.clock() - walltime)

        # Call the finite volume frame construction function from libUtils
        walltime = time.clock()
        self.control_volumes, self.neighbours, self.vor_edges, \
        self.edge_length, maxNgbhs = fvframe.build(lGIDs+1, self.localIDs+1, self.node_coords[:,0], self.node_coords[:,1],
            self.edges[:,:2]+1, self.cells[:,:3]+1, Vor_pts[:,0], Vor_pts[:,1], Vor_edges[:,:2]+1)
        if verbose:
            print(" - construct Finite Volume representation ", time.clock() - walltime)

        # Maximum number of neighbours for each partition
        self.maxNgbh = numpy.array(maxNgbhs)

        return

    def _gather_GIDs(self, lGIDs):
        """
        Gather local IDs to all processors.

        Args:
            lGIDs: numpy integer-type array filled with the global vertex IDs for each local grid located within the partition (including those on the edges).

        Returns
        -------
        exportGIDs
            numpy integer-type array filled with the global vertex IDs ordered by processor ID.
        localPtsNb
            number of points on each local partition.
        """

        # Get global IDs of non-boundary vertex for each TIN
        gid = lGIDs.astype(numpy.int32)
        gids = gid[self.localIDs]
        localPtsNb = len(gids)

        # Find each partition contribution to global dataset
        arraylIDsNb = localPtsNb

        # Gather vertex IDs from each region globally
        exportGIDs = gids

        return exportGIDs, localPtsNb

    def _gather_Area(self, localPtsNb):
        """
        Gather local voronoi area to all processors.

        Args:
            localPtsNb: number of points on each local partition.

        Returns
        -------
        exportVols
            numpy float-type array containing the voronoi area for each TIN node.
        """

        # Get local volume declaration
        vols = numpy.zeros(localPtsNb, dtype=numpy.float)
        vols = self.control_volumes[self.localIDs]
        volsFLT = vols.astype(numpy.float64)

        # Find each partition contribution to global dataset
        arraylocNb = len(volsFLT)

        # Gather flatten neighbour array definition from each region globally
        exportVols = volsFLT

        return exportVols

    def _gather_Neighbours(self, localPtsNb, totPts):
        """
        Gather local neigbours ID to all processors.

        Args:
            localPtsNb: number of points on each local partition.
            totPts: total number of points on the global TIN surface.

        Returns
        -------
        exportNgbhIDs
            numpy integer-type array filled with the global neighbourhood IDs.
        shape
            shape of the neighbours array.
        ngbhNbs
            numpy integer-type array filled with the local neighbourhood IDs.
        """

        # Get local neighbourhood declaration
        ngbh = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.int)
        ngbh.fill(-2)
        ngbh = self.neighbours[self.localIDs,:self.maxNgbh]
        ngbh = numpy.ravel(ngbh)
        ngbhINT = ngbh.astype(numpy.int32)

        # Find each partition contribution to global dataset
        ngbhNbs = len(ngbhINT)

        # Gather flatten neighbour array definition from each region globally
        shape = (totPts, self.maxNgbh)
        #globalNgbh = numpy.zeros(sum(ngbhNbs),dtype=ngbhINT.dtype)
        globalNgbh = ngbhINT
        exportNgbhIDs = numpy.reshape(globalNgbh,shape)

        return exportNgbhIDs, shape, ngbhNbs

    def _gather_Edges(self, localPtsNb, shape, ngbhNbs):
        """
        Gather local edges to all processors.

        Args:
            localPtsNb: number of points on each local partition.
            shape: shape of the neighbours array.
            ngbhNbs: numpy integer-type array filled with the local neighbourhood IDs.

        Returns
        -------
        exportEdges
            numpy float-type array containing the lengths to each neighbour.
        """

        # Get local edges length declaration
        edges = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float)
        edges = self.edge_length[self.localIDs,:self.maxNgbh]
        edges = numpy.ravel(edges)
        edgesFLT = edges.astype(numpy.float64)

        # Gather flatten array definition from each region globally
        #globalEdges = numpy.zeros(sum(ngbhNbs),dtype=edgesFLT.dtype)
        globalEdges = edgesFLT
        exportEdges = numpy.reshape(globalEdges,shape)

        return exportEdges

    def _gather_VorEdges(self, localPtsNb, shape, ngbhNbs):
        """
        Gather local voronoi edges to all processors.

        Args:
            localPtsNb: number of points on each local partition.
            shape: shape of the neighbours array.
            ngbhNbs: numpy integer-type array filled with the local neighbourhood IDs.

        Returns
        -------
        exportVors
            numpy float-type array containing the voronoi edge lengths to each neighbour.
        """

        # Get local voronoi length declaration
        vors = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float)
        vors = self.vor_edges[self.localIDs,:self.maxNgbh]
        vors = numpy.ravel(vors)
        vorsFLT = vors.astype(numpy.float64)

        # Gather flatten array definition from each region globally
        globalVors = vorsFLT
        exportVors = numpy.reshape(globalVors,shape)

        return exportVors

    def construct_FV(self, inIDs, lGIDs, totPts, res, verbose=False):
        """
        Called function to build the Finite Volume discretisation of badlands TIN grid.

        The approach provides an efficient method for storing, accessing, and updating a Delaunay
        triangulation and its associated Voronoi diagram. It is inspired by the method described in
        Tucker et al. (2001).

        Note:
            Tucker et al., 2001: An object-oriented framework for distributed hydrologic and
            geomorphic modeling using triangulated irregular networks, Computers & Geosciences,
            27 (8), 959-973, `doi:10.1016/S0098-3004(00)00134-5`_.

        .. _doi:10.1016/S0098-3004(00)00134-5: https://doi.org/10.1016/S0098-3004(00)00134-5

        Args:
            nIDs: numpy integer-type array filled with the global vertex IDs for each local grid located within the partition (not those on the edges).
            lGIDs: numpy integer-type array filled with the global vertex IDs for each local grid located within the partition (including those on the edges).
            totPts: total number of points on the global TIN surface.
            res: resolution of the tin edges.
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).

        Returns
        -------
        exportGIDs
            numpy integer-type array filled with the global vertex IDs ordered by processor ID.
        exportNgbhIDs
            numpy integer-type array filled with the global neighbourhood IDs.
        exportEdges
            numpy float-type array containing the lengths to each neighbour.
        exportVors
            numpy float-type array containing the voronoi edge lengths to each neighbour.
        exportVols
            numpy float-type array containing the voronoi area for each TIN node.


        Important:
            In the background the code is calling **Triangle** a python wrapper around *Jonathan Richard
            Shewchuk*'s two-dimensional quality mesh generator and delaunay triangulator library,
            available `here <quake_>`_.  The source is available on Github_.

        .. _quake: http://www.cs.cmu.edu/~quake/triangle.html
        .. _Github: https://github.com/drufat/triangle
        """

        # Define each partition nodes global IDs
        inArrays = numpy.in1d(lGIDs, inIDs)
        ids = numpy.where(inArrays == True)[0]
        self.partIDs = numpy.zeros(len(lGIDs))
        self.partIDs.fill(-1)
        self.partIDs[ids] = 0

        # Get each partition local node ID
        self.localIDs = numpy.where(self.partIDs == 0)[0]

        # Call finite volume function
        self._FV_utils(lGIDs)

        # Gather processor dataset together
        walltime = time.clock()

        # Gather vertex IDs from each region globally
        exportGIDs, localPtsNb = self._gather_GIDs(lGIDs)
        # Gather voronoi area from each region globally
        exportVols = self._gather_Area(localPtsNb)
        # Gather neighbourhood IDs from each region globally
        exportNgbhIDs, shape, ngbhNbs = self._gather_Neighbours(localPtsNb, totPts)
        # Gather edges lengths from each region globally
        exportEdges = self._gather_Edges(localPtsNb, shape, ngbhNbs)
        maxdist = numpy.sqrt(2.*res**2)
        exportEdges[exportEdges > 2.*maxdist] = maxdist
        # Gather voronoi edges lengths from each region globally
        exportVors = self._gather_VorEdges(localPtsNb, shape, ngbhNbs)

        return exportGIDs, exportNgbhIDs, exportEdges, exportVors, exportVols
