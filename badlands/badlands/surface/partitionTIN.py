##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Approach to distribute mesh parameters amongst processors.

Warning:
    This set of functions are inherited from a former version of **badlands** when
    a parallel approach was implemented!
"""

import time
import numpy
import tribad as triangle

import os
if 'READTHEDOCS' not in os.environ:
    from badlands import flowalgo

from numpy import random


def simple(X, Y, Xdecomp=1, Ydecomp=1):
    """
    This function defines a simple partitioning of the computational domain based on
    row and column wise decomposition.

    Note:
        The method is relatively fast compared to other techniques
        but lack load-balancing operations.

    The purpose of this function is **twofold**.

    First, it efficiently decomposes the domain based on the number of processors available.
    Second, it returns to all processors the partition IDs for each nodes composing the TIN.

    Args:
        X: numpy array containing the X coordinates of the TIN vertices.
        Y: numpy array containing the Y coordinates of the TIN vertices.
        Xdecomp: integers that specifies the number of processors along X axis (default: 1).
        Ydecomp: integers that specifies the number of processors along Y axis (default: 1).

    Important:
        It is a requirement that the number of processors (:math:`P_{nb}`) used matches the
        proposed decomposition:

        .. math::
           P_{nb} = n_X \cdot n_Y

    Returns
    -------
    partID
        numpy integer-type array filled with the ID of the partition each node belongs to.
    nbprocX
        integer specifying the number of processors along X axis.
    nbprocY
        integer specifying the number of processors along Y axis.
    """

    xmin = X.min()
    xmax = X.max()
    ymin = Y.min()
    ymax = Y.max()
    size = 1
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
    This function defines a simple partitioning of the computational domain based on
    row and column wise decomposition.

    The purpose of the class is to efficiently decomposed the domain from the number of
    processors defined along the X and Y axes. It creates an overlapping region between each
    partition. It also returns to each processor their contained TIN vertices and finally it
    builds a local TIN on each of the decomposed domain

    Args:
        X : numpy arrays containing the X coordinates of the TIN vertices.
        Y : numpy arrays containing the Y coordinates of the TIN vertices.
        nbprocX : integers that specifies the number of processors along X axis.
        nbprocY : integers that specifies the number of processors along Y axis.
        overlapLen : defining the length of the overlapping region.

    Returns
    -------
    globIDs
        numpy integer-type array containing for local nodes their global IDs.
    localTIN
        triangle class defining TIN model.


    Important:
        It is a requirement that the number of processors (:math:`P_{nb}`) used matches the
        proposed decomposition:

        .. math::
           P_{nb} = n_X \cdot n_Y


    Caution:
        The method is relatively fast compared to other techniques but lacks
        *load-balancing* operations.

    """

    walltime = time.clock()
    size = 1
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
    partID = flowalgo.overlap(X,Y,Xst[0],Yst[0],Xed[0],Yed[0])

    # Extract local domain nodes global ID
    globIDs = numpy.where(partID > -1)[0]

    # Build local TIN
    data = numpy.column_stack((X,Y))
    localTIN = triangle.triangulate( {'vertices':data[globIDs,:2]},' ')

    if verbose:
        print(" - partition TIN including shadow zones ", time.clock() - walltime)

    return globIDs, localTIN
