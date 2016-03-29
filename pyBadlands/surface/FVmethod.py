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
    This class builds parameters required for the Finite Volume mesh algorithm.
    
    It creates the following for each node: 
        1. the voronoi cell area
        2. an ordered list of voronoi edges lenght
        
    Parameters
    ----------
    string : nodes
        The 2D coordinates of TIN nodes coordinates.
        
    string : cells
        IDs of each node defining TIN's cells.
        
    variable: edges
        IDs of each edges from the TIN.
    """
    
    def __init__(self, nodes, cells, edges):

        self.node_coords = nodes
        self.edges = edges
        self.cells = cells
        self.control_volumes = None
        self.neighbours = None
        self.vor_edges = None
        self.edge_lenght = None
        self.fillH = None
        self.partIDs = None
        self.maxNgbh = None
        self.localIDs = None
        
        return
    
    def _FV_utils(self, allIDs):
        """
        This function constructs the Finite Volume discretisation for each local
        triangularised grid.
        
        Parameters
        ----------
        variable : allIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).
            
        variable: inIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (not those on the edges).
        """
        
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        
        # Build the voronoi diagram (dual of the delaunay)
        walltime = time.clock()
        Vor_pts, Vor_edges = triangle.voronoi(self.node_coords)
        if rank == 0:
            print " - build the voronoi diagram ", time.clock() - walltime
        
        # Call the finite volume frame construction function from libUtils     
        walltime = time.clock()
        self.control_volumes, self.neighbours, self.vor_edges, \
        self.edge_lenght, maxNgbhs = FVframe.discretisation.build( \
             allIDs+1, self.localIDs+1, self.node_coords[:,0], self.node_coords[:,1],
            self.edges[:,:2]+1, self.cells[:,:3]+1, Vor_pts[:,0], Vor_pts[:,1], \
            Vor_edges[:,:2]+1)
        if rank == 0:
            print " - construct Finite Volume representation ", time.clock() - walltime
        
        # Maximum number of neighbours for each partition
        maxNgbh = numpy.array(maxNgbhs)
        comm.Allreduce(mpi.IN_PLACE,maxNgbh,op=mpi.MAX)
        self.maxNgbh = maxNgbh
        
        return 
    
    def _gather_GIDs(self, allIDs):
        """
        Gather local IDs to all processors.
        
        Parameters
        ----------
        variable : allIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).
            
        Return
        ----------
        variable: exportGIDs
            Numpy integer-type array filled with the global vertex IDs ordered by processor ID.

        variable : localPtsNb
            Number of points on each local partition.
        """
        
        comm = mpi.COMM_WORLD
        
        # Get global IDs of non-boundary vertex for each TIN
        gid = allIDs.astype(numpy.int32)
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
        variable : localPtsNb
            Number of points on each local partition.
            
        Return
        ----------
        variable : exportVols
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
        variable : localPtsNb
            Number of points on each local partition.
            
        variable : totPts
            Total number of points on the global TIN surface.
            
        Return
        ----------
        variable : exportNgbhIDs
            Numpy integer-type array filled with the global neighbourhood IDs.
            
        variable : shape
            Shape of the neighbours array.
            
        variable : ngbhNbs
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
        variable : localPtsNb
            Number of points on each local partition.
            
        variable : shape
            Shape of the neighbours array.
            
        variable : ngbhNbs
            Numpy integer-type array filled with the local neighbourhood IDs.
            
        Return
        ----------
        variable : exportEdges
            Numpy float-type array containing the lenghts to each neighbour.
        """
        
        comm = mpi.COMM_WORLD
        
        # Get local edges length declaration
        edges = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float) 
        edges = self.edge_lenght[self.localIDs,:self.maxNgbh]
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
        variable : localPtsNb
            Number of points on each local partition.
            
        variable : shape
            Shape of the neighbours array.
            
        variable : ngbhNbs
            Numpy integer-type array filled with the local neighbourhood IDs.
            
        Return
        ----------
        variable : exportVors
            Numpy float-type array containing the voronoi edge lenghts to each neighbour.    
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
    
    def construct_FV(self, inIDs, allIDs, totPts):
        """
        Called function to build the Finite Volume discretisation of Badlands TIN grid.
        
        Parameters
        ----------
        variable: inIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (not those on the edges).

        variable : allIDs
            Numpy integer-type array filled with the global vertex IDs for each local grid located
            within the partition (including those on the edges).
            
        variable : totPts
            Total number of points on the global TIN surface.
            
        Return
        ----------
        variable: exportGIDs
            Numpy integer-type array filled with the global vertex IDs ordered by processor ID.

        variable : exportNgbhIDs
            Numpy integer-type array filled with the global neighbourhood IDs.
            
        variable : exportEdges
            Numpy float-type array containing the lenghts to each neighbour.
            
        variable : exportVors
            Numpy float-type array containing the voronoi edge lenghts to each neighbour.    
            
        variable : exportVols
            Numpy float-type array containing the voronoi area for each TIN node.    
        """
        
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        
        # Define each partition nodes global IDs
        inArrays = numpy.in1d(allIDs, inIDs)
        ids = numpy.where(inArrays == True)[0]
        self.partIDs = numpy.zeros(len(allIDs)) 
        self.partIDs.fill(-1)
        self.partIDs[ids] = rank
        
        # Get each partition local node ID
        self.localIDs = numpy.where(self.partIDs == rank)[0]
        
        # Call finite volume function
        self._FV_utils(allIDs)
        
        # Gather processor dataset together
        walltime = time.clock()
        
        # Gather vertex IDs from each region globally
        exportGIDs, localPtsNb = self._gather_GIDs(allIDs)        
        # Gather voronoi area from each region globally
        exportVols = self._gather_Area(localPtsNb)     
        # Gather neighbourhood IDs from each region globally
        exportNgbhIDs, shape, ngbhNbs = self._gather_Neighbours(localPtsNb, totPts)
        # Gather edges lenghts from each region globally
        exportEdges = self._gather_Edges(localPtsNb, shape, ngbhNbs)
        # Gather voronoi edges lenghts from each region globally
        exportVors = self._gather_VorEdges(localPtsNb, shape, ngbhNbs)
        
        if rank == 0:
            print " - perform MPI communication ", time.clock() - walltime
            
        # Get global IDs of non-boundary vertex for each TIN
        #gid = allIDs.astype(numpy.int32)
        #gids = gid[localIDs]
        #localPtsNb = len(gids)
        
        # Find each partition contribution to global dataset
        #arraylIDsNb = comm.allgather(localPtsNb)

        # Gather vertex IDs from each region globally
        #exportGIDs = numpy.zeros(sum(arraylIDsNb),dtype=gids.dtype)
        #comm.Allgatherv(sendbuf=[gids, mpi.INTEGER4], 
        #    recvbuf=[exportGIDs, (arraylIDsNb, None), mpi.INTEGER4]) 

        
        # Get local neighbourhood declaration
        #ngbh = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.int) 
        #ngbh.fill(-2)
        #ngbh = self.neighbours[localIDs,:self.maxNgbh]
        #ngbh = numpy.ravel(ngbh)
        #ngbhINT = ngbh.astype(numpy.int32)

        # Find each partition contribution to global dataset
        #ngbhNbs = comm.allgather(len(ngbhINT))

        # Gather flatten neighbour array definition from each region globally
        #shape1 = (totPts, self.maxNgbh)
        #globalNgbh = numpy.zeros(sum(ngbhNbs),dtype=ngbhINT.dtype)
        #comm.Allgatherv(sendbuf=[ngbhINT, mpi.INT], 
        #             recvbuf=[globalNgbh, (ngbhNbs, None), mpi.INT]) 
        #exportNgbhIDs = numpy.reshape(globalNgbh,shape1)
        
        
        ################
        
        # Get local edges length declaration
        #edges = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float) 
        #edges = self.edge_lenght[localIDs,:self.maxNgbh]
        #edges = numpy.ravel(edges)
        #edgesFLT = edges.astype(numpy.float32)
        
        # Gather flatten array definition from each region globally
        #globalEdges = numpy.zeros(sum(ngbhNbs),dtype=edgesFLT.dtype)
        #comm.Allgatherv(sendbuf=[edgesFLT, mpi.FLOAT], 
        #             recvbuf=[globalEdges, (ngbhNbs, None), mpi.FLOAT]) 
        #exportEdges = numpy.reshape(globalEdges,shape1)
        
        
        ################
        
        # Get local voronoi length declaration
        #vors = numpy.zeros((localPtsNb,self.maxNgbh), dtype=numpy.float) 
        #vors = self.vor_edges[localIDs,:self.maxNgbh]
        #vors = numpy.ravel(vors)
        #vorsFLT = vors.astype(numpy.float32)
        
        # Gather flatten array definition from each region globally
        #globalVors = numpy.zeros(sum(ngbhNbs),dtype=vorsFLT.dtype)
        #comm.Allgatherv(sendbuf=[vorsFLT, mpi.FLOAT], 
        #             recvbuf=[globalVors, (ngbhNbs, None), mpi.FLOAT]) 
        #exportVors = numpy.reshape(globalVors,shape1)
        
        # Get local volume declaration
        #vols = numpy.zeros(localPtsNb, dtype=numpy.float) 
        #vols = self.control_volumes[localIDs]
        #volsFLT = vols.astype(numpy.float32)
        
        # Find each partition contribution to global dataset
        #ngbhNbs2 = comm.allgather(len(volsFLT))

        # Gather flatten neighbour array definition from each region globally
        #exportVols = numpy.zeros(sum(ngbhNbs2),dtype=volsFLT.dtype)
        #comm.Allgatherv(sendbuf=[volsFLT, mpi.FLOAT], 
        #             recvbuf=[exportVols, (ngbhNbs2, None), mpi.FLOAT])  
        
        #if rank == 0:
        #    print " - perform MPI communication ", time.clock() - walltime
        
        return exportGIDs, exportNgbhIDs, exportEdges, exportVors, exportVols

