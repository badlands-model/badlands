##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to Badlands SP finite volume discretisation.
"""

import time
import numpy
from pyBadlands.libUtils import FVframe
import warnings
import triangle
import mpi4py.MPI as mpi


class FVmethod:
    """
    This class builds paramters required for the Finite Volume mesh algorithm.

    It creates the following for each node:
        1. the voronoi cell area
        2. an ordered list of voronoi edges length

    Parameters
    ----------
    nodes : string
        The 2D coordinates of TIN nodes coordinates.

    cells : string
        IDs of each node defining TIN's cells.

    edges
        IDs of each edges from the TIN.
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
        This function constructs the Finite Volume discretisation for each local
        triangularised grid.

        Parameters
        ----------
        lGIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).

        inIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (not those on the edges).
        """

        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()

        # Build the voronoi diagram (dual of the delaunay)
        walltime = time.clock()
        Vor_pts, Vor_edges = triangle.voronoi(self.node_coords)
        if rank == 0 and verbose:
            print " - build the voronoi diagram ", time.clock() - walltime

        # Call the finite volume frame construction function from libUtils
        walltime = time.clock()
        self.control_volumes, self.neighbours, self.vor_edges, \
        self.edge_length, maxNgbhs = FVframe.discretisation.build( \
             lGIDs+1, self.localIDs+1, self.node_coords[:,0], self.node_coords[:,1],
            self.edges[:,:2]+1, self.cells[:,:3]+1, Vor_pts[:,0], Vor_pts[:,1], \
            Vor_edges[:,:2]+1)
        if rank == 0 and verbose:
            print " - construct Finite Volume representation ", time.clock() - walltime

        # Maximum number of neighbours for each partition
        maxNgbh = numpy.array(maxNgbhs)
        comm.Allreduce(mpi.IN_PLACE,maxNgbh,op=mpi.MAX)
        self.maxNgbh = maxNgbh

    def _gather_GIDs(self, lGIDs):
        """
        Gather local IDs to all processors.

        Parameters
        ----------
        lGIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).

        Returns
        -------
        exportGIDs
            Numpy integer-type array filled with the global vertex IDs ordered by processor ID.

        localPtsNb
            Number of points on each local partition.
        """

        comm = mpi.COMM_WORLD

        # Get global IDs of non-boundary vertex for each TIN
        gid = lGIDs.astype(numpy.int32)
        gids = gid[self.localIDs]
        localPtsNb = len(gids)

        # Find each partition contribution to global dataset
        arraylIDsNb = comm.allgather(localPtsNb)

        # Gather vertex IDs from each region globally
        exportGIDs = numpy.zeros(sum(arraylIDsNb),dtype=gids.dtype)
        comm.Allgatherv(sendbuf=[gids, mpi.INTEGER4],
            recvbuf=[exportGIDs, (arraylIDsNb, None), mpi.INTEGER4])

        return exportGIDs, localPtsNb

    def _gather_Area(self, localPtsNb):
        """
        Gather local voronoi area to all processors.

        Parameters
        ----------
        localPtsNb
            Number of points on each local partition.

        Returns
        -------
        exportVols
            Numpy float-type array containing the voronoi area for each TIN node.
        """

        comm = mpi.COMM_WORLD

        # Get local volume declaration
        vols = numpy.zeros(localPtsNb, dtype=numpy.float)
        vols = self.control_volumes[self.localIDs]
        volsFLT = vols.astype(numpy.float32)

        # Find each partition contribution to global dataset
        arraylocNb = comm.allgather(len(volsFLT))

        # Gather flatten neighbour array definition from each region globally
        exportVols = numpy.zeros(sum(arraylocNb),dtype=volsFLT.dtype)
        comm.Allgatherv(sendbuf=[volsFLT, mpi.FLOAT],
                     recvbuf=[exportVols, (arraylocNb, None), mpi.FLOAT])

        return exportVols

    def _gather_Neighbours(self, localPtsNb, totPts):
        """
        Gather local neigbours ID to all processors.

        Parameters
        ----------
        localPtsNb
            Number of points on each local partition.

        totPts
            Total number of points on the global TIN surface.

        Returns
        -------
        exportNgbhIDs
            Numpy integer-type array filled with the global neighbourhood IDs.

        shape
            Shape of the neighbours array.

        ngbhNbs
            Numpy integer-type array filled with the local neighbourhood IDs.
        """

        comm = mpi.COMM_WORLD

        # Get local neighbourhood declaration
        ngbh = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.int)
        ngbh.fill(-2)
        ngbh = self.neighbours[self.localIDs,:self.maxNgbh]
        ngbh = numpy.ravel(ngbh)
        ngbhINT = ngbh.astype(numpy.int32)

        # Find each partition contribution to global dataset
        ngbhNbs = comm.allgather(len(ngbhINT))

        # Gather flatten neighbour array definition from each region globally
        shape = (totPts, self.maxNgbh)
        globalNgbh = numpy.zeros(sum(ngbhNbs),dtype=ngbhINT.dtype)
        comm.Allgatherv(sendbuf=[ngbhINT, mpi.INT],
                     recvbuf=[globalNgbh, (ngbhNbs, None), mpi.INT])
        exportNgbhIDs = numpy.reshape(globalNgbh,shape)

        return exportNgbhIDs, shape, ngbhNbs

    def _gather_Edges(self, localPtsNb, shape, ngbhNbs):
        """
        Gather local edges to all processors.

        Parameters
        ----------
        localPtsNb
            Number of points on each local partition.

        shape
            Shape of the neighbours array.

        ngbhNbs
            Numpy integer-type array filled with the local neighbourhood IDs.

        Returns
        -------
        exportEdges
            Numpy float-type array containing the lengths to each neighbour.
        """

        comm = mpi.COMM_WORLD

        # Get local edges length declaration
        edges = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float)
        edges = self.edge_length[self.localIDs,:self.maxNgbh]
        edges = numpy.ravel(edges)
        edgesFLT = edges.astype(numpy.float32)

        # Gather flatten array definition from each region globally
        globalEdges = numpy.zeros(sum(ngbhNbs),dtype=edgesFLT.dtype)
        comm.Allgatherv(sendbuf=[edgesFLT, mpi.FLOAT],
                     recvbuf=[globalEdges, (ngbhNbs, None), mpi.FLOAT])
        exportEdges = numpy.reshape(globalEdges,shape)

        return exportEdges

    def _gather_VorEdges(self, localPtsNb, shape, ngbhNbs):
        """
        Gather local voronoi edges to all processors.

        Parameters
        ----------
        localPtsNb
            Number of points on each local partition.

        shape
            Shape of the neighbours array.

        ngbhNbs
            Numpy integer-type array filled with the local neighbourhood IDs.

        Returns
        -------
        exportVors
            Numpy float-type array containing the voronoi edge lengths to each neighbour.
        """

        comm = mpi.COMM_WORLD

        # Get local voronoi length declaration
        vors = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float)
        vors = self.vor_edges[self.localIDs,:self.maxNgbh]
        vors = numpy.ravel(vors)
        vorsFLT = vors.astype(numpy.float32)

        # Gather flatten array definition from each region globally
        globalVors = numpy.zeros(sum(ngbhNbs),dtype=vorsFLT.dtype)
        comm.Allgatherv(sendbuf=[vorsFLT, mpi.FLOAT],
                     recvbuf=[globalVors, (ngbhNbs, None), mpi.FLOAT])
        exportVors = numpy.reshape(globalVors,shape)

        return exportVors

    def construct_FV(self, inIDs, lGIDs, totPts, res, verbose=False):
        """
        Called function to build the Finite Volume discretisation of Badlands TIN grid.

        Parameters
        ----------
        nIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (not those on the edges).

        lGIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).

        totPts
            Total number of points on the global TIN surface.

        res
            Resolution of the tin edges.

        Returns
        -------
        xportGIDs
            Numpy integer-type array filled with the global vertex IDs ordered by processor ID.

        exportNgbhIDs
            Numpy integer-type array filled with the global neighbourhood IDs.

        exportEdges
            Numpy float-type array containing the lengths to each neighbour.

        exportVors
            Numpy float-type array containing the voronoi edge lengths to each neighbour.

        exportVols
            Numpy float-type array containing the voronoi area for each TIN node.
        """

        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()

        # Define each partition nodes global IDs
        inArrays = numpy.in1d(lGIDs, inIDs)
        ids = numpy.where(inArrays == True)[0]
        self.partIDs = numpy.zeros(len(lGIDs))
        self.partIDs.fill(-1)
        self.partIDs[ids] = rank

        # Get each partition local node ID
        self.localIDs = numpy.where(self.partIDs == rank)[0]

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

        if rank == 0 and verbose:
            print " - perform MPI communication ", time.clock() - walltime

        return exportGIDs, exportNgbhIDs, exportEdges, exportVors, exportVols
