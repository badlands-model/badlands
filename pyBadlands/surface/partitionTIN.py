##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module proposes 2 methods of triangular irregular network (TIN) partitioning.
"""

import time
import numpy
import triangle
from pyBadlands.libUtils import FASTloop
import mpi4py.MPI as mpi

from numpy import random

def get_closest_factors(size):
    """
    This function finds the two closest integers which, when multiplied, equal a given number.
    This is used to defined the partition of the regular and TIN grids.

    Parameters
    ----------
    size
        Integer corresponding to the number of CPUs that are used.

    Returns
    -------
    partID
        Numpy integer-type array filled with the ID of the partition each node belongs to.

    nb1
        Integers which specify the number of processors along each axis.

    nb2
        Integers which specify the number of processors along each axis.
    """
    factors =  []
    for i in range(1, size + 1):
        if size % i == 0:
            factors.append(i)
    factors = numpy.array(factors)
    if len(factors)%2 == 0:
        n1 = int( len(factors)/2 ) - 1
        n2 = n1 + 1
    else:
        n1 = int( len(factors)/2 )
        n2 = n1
    nb1 = factors[n1]
    nb2 = factors[n2]
    if nb1*nb2 != size:
        raise ValueError('Error in the decomposition grid: the number of domains \
        decomposition is not equal to the number of CPUs allocated')

    return nb1, nb2

def simple(X, Y, Xdecomp=1, Ydecomp=1):
    """
    This function defines a simple partitioning of the computational domain based on
    row and column wise decomposition. The method is relatively fast compared to other techniques
    but lack load-balancing operations.

    The purpose of the class is:

        1. to efficiently decomposed the domain from the number of processors defined along the X and Y axes
        2. to return to all processors the partition IDs for each vertice of the TIN

    Parameters
    ----------
    X
        Numpy array containing the X coordinates of the TIN vertices.

    Y
        Numpy array containing the Y coordinates of the TIN vertices.

    Xdecomp
        Integers which specify the number of processors along each axis. It is a requirement that
        the number of processors used matches the proposed decomposition:

            >> nb CPUs = nbprocX x nbprocY

    Ydecomp
        Integers which specify the number of processors along each axis. It is a requirement that
        the number of processors used matches the proposed decomposition:

            >> nb CPUs = nbprocX x nbprocY

    Returns
    -------
    partID
        Numpy integer-type array filled with the ID of the partition each node belongs to.

    nbprocX
        Integers which specify the number of processors along each axis.

    nbprocY
        Integers which specify the number of processors along each axis.
    """

    # Initialise MPI communications
    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    xmin = X.min()
    xmax = X.max()
    ymin = Y.min()
    ymax = Y.max()

    if Xdecomp == 1 and Ydecomp == 1 and size > 1:
        n1,n2 = get_closest_factors(size)
        if xmax-xmin > ymax-ymin :
            nbprocX = n2
            nbprocY = n1
        else:
            nbprocX = n1
            nbprocY = n2
    else:
        nbprocX = Xdecomp
        nbprocY = Ydecomp

    # Check decomposition versus CPUs number
    if size != nbprocX*nbprocY:
        raise ValueError('Error in the decomposition grid: the number of domains \
        decomposition is not equal to the number of CPUs allocated')

    # Define output type and size
    partID = numpy.zeros( len(X), dtype=numpy.uint32 )

    # Get extent of X partition
    nbX = int((xmax-xmin)/nbprocX)
    Xstart = numpy.zeros( nbprocX )
    Xend = numpy.zeros( nbprocX )
    for p in range(nbprocX):
        Xstart[p] = p*nbX+xmin
        Xend[p] = Xstart[p]+nbX
    Xend[nbprocX-1]=xmax

    # Get extent of Y partition
    nbY = int((ymax-ymin)/nbprocY)
    Ystart = numpy.zeros( nbprocY )
    Yend = numpy.zeros( nbprocY )
    for p in range(nbprocY):
        Ystart[p] = p*nbY+ymin
        Yend[p] = Ystart[p]+nbY
    Yend[nbprocY-1]=ymax

    # Fill partition ID based on node coordinates
    for id in range(len(X)):
        ix = 0
        pX = Xend[ix]
        while X[id] > pX :
            ix += 1
            pX = Xend[ix]
        iy = 0
        pY = Yend[iy]
        while Y[id] > pY :
            iy += 1
            pY = Yend[iy]
        partID[id] = ix + iy * nbprocX

    return partID, nbprocX, nbprocY

def overlap(X, Y, nbprocX, nbprocY, overlapLen, verbose=False):
    """
    This function defines a simple partitioning of the computational domain based on row and column
    wise decomposition and add an overlap between each domain. The method is relatively fast compared
    to other techniques but lack load-balancing operations.

    The purpose of the class is:

        1. to efficiently decomposed the domain from the number of processors defined along the X and Y axes
        2. to create an overlapping region between each partition
        3. to return to each processor their contained TIN vertices
        4. to built a local TIN on each of the decomposed domain

    Parameters
    ----------
    X
        Numpy arrays containing the X coordinates of the TIN vertices.

    Y
        Numpy arrays containing the Y coordinates of the TIN vertices.

    nbprocX
        Integers which specify the number of processors along each axis. It is a requirement that
        the number of processors used matches the proposed decomposition:

            >> nb CPUs = nbprocX x nbprocY

    nbprocY
        Integers which specify the number of processors along each axis. It is a requirement that
        the number of processors used matches the proposed decomposition:

            >> nb CPUs = nbprocX x nbprocY

    overlapLen
        Float defining the length of the overlapping region.

    Returns
    -------
    globIDs
        Numpy integer-type array containing for local nodes their global IDs.

    localTIN
        Triangle class representing local TIN coordinates and parameters.
    """

    # Initialise MPI communications
    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    walltime = time.clock()

    # Check decomposition versus CPUs number
    if size != nbprocX*nbprocY:
        raise ValueError('Error in the decomposition grid: the number of domains \
        decomposition is not equal to the number of CPUs allocated')

    # Get extent of X partition
    xmin = X.min()
    xmax = X.max()
    nbX = int((xmax-xmin)/nbprocX)
    Xstart = numpy.zeros( nbprocX )
    Xend = numpy.zeros( nbprocX )
    for p in range(nbprocX):
        if p == 0:
            Xstart[p] = p*nbX+xmin
            Xend[p] = Xstart[p]+nbX+overlapLen
        else:
            Xstart[p] = p*nbX+xmin-overlapLen
            Xend[p] = Xstart[p]+nbX+2*overlapLen
    Xend[nbprocX-1]=xmax

    # Get extent of Y partition
    ymin = Y.min()
    ymax = Y.max()
    nbY = int((ymax-ymin)/nbprocY)
    Ystart = numpy.zeros( nbprocY )
    Yend = numpy.zeros( nbprocY )
    for p in range(nbprocY):
        if p == 0:
            Ystart[p] = p*nbY+ymin
            Yend[p] = Ystart[p]+nbY+overlapLen
        else:
            Ystart[p] = p*nbY+ymin-overlapLen
            Yend[p] = Ystart[p]+nbY+2*overlapLen
    Yend[nbprocY-1]=ymax

    # Define partitions ID globally
    Xst = numpy.zeros( size )
    Xed = numpy.zeros( size )
    Yst = numpy.zeros( size )
    Yed = numpy.zeros( size )
    for q in range(nbprocX):
        for p in range(nbprocY):
            Xst[q + p * nbprocX] = Xstart[q]
            Yst[q + p * nbprocX] = Ystart[p]
            Xed[q + p * nbprocX] = Xend[q]
            Yed[q + p * nbprocX] = Yend[p]

    # Loop over node coordinates and find if they belong to local partition
    # Note: used a Cython/Fython class to increase search loop performance... in libUtils
    partID = FASTloop.part.overlap(X,Y,Xst[rank],Yst[rank],Xed[rank],Yed[rank])

    # Extract local domain nodes global ID
    globIDs = numpy.where(partID > -1)[0]

    # Build local TIN
    data = numpy.column_stack((X,Y))
    localTIN = triangle.triangulate(dict(vertices=data[globIDs,:2]),' ')

    if rank == 0 and verbose:
        print " - partition TIN including shadow zones ", time.clock() - walltime

    return globIDs, localTIN

def _robin_distribution(X,Y):
    """
    This function defines an initial distribution using round-robin algorithm.

    Parameters
    ----------
    X: Numpy array containing the X coordinates of the TIN vertices.

    Y: Numpy array containing the Y coordinates of the TIN vertices.

    Returns
    -------
    GIDs
        Numpy integer-type array containing local nodes global IDs.

    Lx
        Numpy float-type array containing local nodes X coordinates.

    Ly
        Numpy float-type array containing local nodes Y coordinates.
    """

    # Initialise MPI communications
    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    numGlob = len(X)
    numLoc = 0

    # Round robin initial distribution
    for i in xrange(numGlob):
        if i%size == rank:
            numLoc += 1

    # Allocate data for round robin distribution grid
    GIDs = numpy.zeros(numLoc,dtype=int)
    Lx = numpy.zeros(numLoc,dtype=float)
    Ly = numpy.zeros(numLoc,dtype=float)

    # Fill data for round robin distribution grid
    idx = 0
    for i in xrange(numGlob):
        # Assumes gids start at 1, gives round robin initial distribution
        if i%size == rank:
            GIDs[idx] = i
            Lx[idx] = X[i]
            Ly[idx] = Y[i]
            idx += 1

    return GIDs, Lx, Ly

def _compute_partition_ghosts(size, neighbours, partID):
    """
    This function find the ghosts (nodes) in the vicinity of each decomposition zone.

    Parameters
    ----------
    size
        Number of processors.

    neighbours
        Numpy integer-type array containing for each nodes its neigbhours IDs

    partID
        Numpy integer-type array containing for each nodes its partition ID

    Returns
    -------
    ghosts
        List containing the ghost nodes for each partition.
    """

    ghosts = {}
    for p in range(size):
        ghostIDs = numpy.array([], dtype=int)
        localIDs = numpy.where(partID == p)[0]
        for id in range(len(localIDs)):
            ids = numpy.where( (partID[neighbours[localIDs[id]]] != p) & (neighbours>=0))[0]
            for k in range(len(ids)):
                ghostIDs = numpy.append(ghostIDs, neighbours[localIDs[id]][ids[k]])
        ghosts[p] = numpy.unique(ghostIDs)

    return ghosts
