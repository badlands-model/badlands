##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file contains functions to evaluate elevation on the *computational irregular grid*.

.. image:: img/tin.png
   :scale: 70 %
   :alt: TIN grid
   :align: center

From the initial regular grid, a triangular irregular network (**TIN**) is automatically created
with a resolution that is as high as the regular grid one but that could be down-sampled if the
<resfactor> is used in the input XML file under the :code:`grid structure` module.

.. note::
    Look at the full documentation for additional information.
"""

import time
import numpy

import os
if 'READTHEDOCS' not in os.environ:
    from badlands import pdalgo

from scipy.interpolate import interpn
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator

def _boundary_elevation(elevation, neighbours, edge_length, boundPts, btype):
    """
    This function defines the elevation of the TIN surface edges for 2 different types of conditions:

        * Infinitely flat condition,
        * Continuous slope condition.

    Args:
        elevation: Numpy arrays containing the internal nodes elevation.
        neighbours: Numpy integer-type array containing for each nodes its neigbhours IDs.
        edge_length: Numpy float-type array containing the lengths to each neighbour.
        boundPts: Number of nodes on the edges of the TIN surface.
        btype: Integer defining the type of boundary: 0 for flat and 1 for slope condition.

    Returns:
        - elevation - numpy array containing the updated elevations on the edges.
    """

    # Flat/fixed/wall
    if btype == 0:
        missedPts = []
        for id in range(boundPts):
            ngbhs = neighbours[id,:]
            ids = numpy.where(ngbhs >= boundPts)[0]
            if len(ids) == 1:
                elevation[id] = elevation[ngbhs[ids]]
            elif len(ids) > 1:
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                elevation[id] = elevation[ngbhs[ids[picked]]]
            else:
                missedPts = numpy.append(missedPts,id)

        if len(missedPts) > 0 :
            for p in range(len(missedPts)):
                id = int(missedPts[p])
                ngbhs = neighbours[id,:]
                ids = numpy.where((elevation[ngbhs] < 9.e6) & (ngbhs >= 0))[0]
                if len(ids) == 0:
                    raise ValueError('Error while getting boundary elevation for point ''%d''.' % id)
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                elevation[id] = elevation[ngbhs[ids[picked]]]

    # Slope
    elif btype == 1:
        missedPts = []
        for id in range(boundPts):
            ngbhs = neighbours[id,:]
            ids = numpy.where(ngbhs >= boundPts)[0]
            if len(ids) == 1:
                # Pick closest non-boundary vertice
                ln1 = edge_length[id,ids[0]]
                id1 = ngbhs[ids[0]]
                # Pick closest non-boundary vertice to first picked
                ngbhs2 = neighbours[id1,:]
                ids2 = numpy.where(ngbhs2 >= boundPts)[0]
                lselect = edge_length[id1,ids2]
                if len(lselect) > 0:
                    picked = numpy.argmin(lselect)
                    id2 = ngbhs2[ids2[picked]]
                    ln2 = lselect[picked]
                    elevation[id] = (elevation[id1]-elevation[id2])*(ln2+ln1)/ln2 + elevation[id2]
                else:
                    missedPts = numpy.append(missedPts,id)
            elif len(ids) > 1:
                # Pick closest non-boundary vertice
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                id1 = ngbhs[ids[picked]]
                ln1 = lselect[picked]
                # Pick closest non-boundary vertice to first picked
                ngbhs2 = neighbours[id1,:]
                ids2 = numpy.where(ngbhs2 >= boundPts)[0]
                lselect2 = edge_length[id1,ids2]
                if len(lselect2) > 0:
                    picked2 = numpy.argmin(lselect2)
                    id2 = ngbhs2[ids2[picked2]]
                    ln2 = lselect2[picked2]
                    elevation[id] = (elevation[id1]-elevation[id2])*(ln2+ln1)/ln2 + elevation[id2]
                else:
                    missedPts = numpy.append(missedPts,id)
            else:
                missedPts = numpy.append(missedPts,id)

        if len(missedPts) > 0 :
            for p in range(0,len(missedPts)):
                id = int(missedPts[p])
                ngbhs = neighbours[id,:]
                ids = numpy.where((elevation[ngbhs] < 9.e6) & (ngbhs >= 0))[0]
                if len(ids) == 0:
                    raise ValueError('Error while getting boundary elevation for point ''%d''.' % id)
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                elevation[id] = elevation[ngbhs[ids[picked]]]
        elevation[:boundPts] -= 0.5

    # Associate TIN edge point to the border for ero/dep updates
    parentID = numpy.zeros(boundPts,dtype=int)
    missedPts = []
    for id in range(boundPts):
        ngbhs = neighbours[id,:]
        ids = numpy.where(ngbhs >= boundPts)[0]
        if len(ids) == 1:
            parentID[id] = ngbhs[ids]
        elif len(ids) > 1:
            lselect = edge_length[id,ids]
            picked = numpy.argmin(lselect)
            parentID[id] = ngbhs[ids[picked]]
        else:
            missedPts = numpy.append(missedPts,id)

    if len(missedPts) > 0 :
        for p in range(len(missedPts)):
            id = int(missedPts[p])
            ngbhs = neighbours[id,:]
            ids = numpy.where((elevation[ngbhs] < 9.e6) & (ngbhs >= 0))[0]
            if len(ids) == 0:
                raise ValueError('Error while getting boundary elevation for point ''%d''.' % id)
            lselect = edge_length[id,ids]
            picked = numpy.argmin(lselect)
            parentID[id] = ngbhs[ids[picked]]

    return elevation, parentID

def update_border_elevation(elev, neighbours, edge_length, boundPts, btype='flat'):
    """
    This function computes the boundary elevation based on 3 different conditions:

    1. infinitely flat,
    2. continuous slope ,
    3. wall boundary (closed domain).

    Note:
        The boundary condition are defined via the input XML file.

    Args:
        elev: numpy arrays containing the internal nodes elevation.
        neighbours: numpy integer-type array containing for each nodes its neigbhours IDs.
        edge_length: numpy float-type array containing the lengths to each neighbour.
        boundPts: number of nodes on the edges of the TIN surface.
        btype: integer defining the type of boundary (default: 'flat').

    Returns
    -------
    newelev
        numpy array containing the updated elevations on the edges.
    parentID
        numpy array containing the indices of the associated *inside* node to each boundary node.
    """

    newelev = elev

    if btype == 'wall' or btype == 'flat' or btype == 'slope' or btype == 'fixed' or btype == 'outlet' or btype == 'wall1':
        newelev[:boundPts] = 1.e7
        thetype = 0
        if btype == 'slope' or btype == 'outlet' or btype == 'wall1':
            thetype = 1
        newelev, parentID = _boundary_elevation(elev, neighbours, edge_length, boundPts, thetype)
        if btype == 'wall':
            newelev[:boundPts] = 1.e7
        if btype == 'outlet':
            newelev[1:boundPts] = 1.e7
    else:
        raise ValueError('Unknown boundary type ''%s''.' % btype)

    return newelev, parentID

def getElevation(rX, rY, rZ, coords, interp='linear'):
    """
    This function interpolates elevation from the regular grid to the triamgular mesh
    using **SciPy** *interpn* funciton.

    Args:
        rX: numpy arrays containing the X coordinates from the regular grid.
        rY: numpy arrays containing the Y coordinates from the regular grid.
        rZ: numpy arrays containing the Z coordinates from the regular grid.
        coords: numpy float-type array containing X, Y coordinates for the TIN nodes.
        interp: interpolation method as in *SciPy interpn function* (default: 'linear')

    Returns:
        - elev - numpy array containing the updated elevations for the local domain.
    """

    # Set new elevation to 0
    elev = numpy.zeros(len(coords[:,0]))

    # Get the TIN points elevation values using the regular grid dataset
    elev = interpn( (rX, rY), rZ, (coords[:,:2]), method=interp)

    return elev

def assign_parameter_pit(neighbours, area, diffnb, prop, propa, propb, boundPts, fillTH=1., epsilon=1.e-6):
    """
    This function defines the global variables used in the **pit filling algorithm** described in
    the :code:`pit_stack` function_.

    Args:
        neighbours: numpy integer-type array containing for each nodes its neigbhours IDs.
        edge_length: numpy float-type array containing the length to each neighbour.
        area: numpy float-type array containing the area of each cell.
        diffnb: marine diffusion distribution steps.
        prop: proportion of marine sediment deposited on downstream nodes.
        boundPts: number of nodes on the edges of the TIN surface.
        fillTH: limit the filling algorithm to a specific height to prevent complete filling of depression (default: 1. m).
        epsilon: force the creation of a minimal slope instead of a flat area to ensure continuous flow pathways (default: 1.e-6 m).

    .. _function: https://badlands.readthedocs.io/en/latest/api.html#surface.elevationTIN.pit_stack

    Tip:
        Planchon and Darboux, 2001: A Fast, Simple and Versatile Algorithm to Fill the Depressions of Digital Elevation Models - Catena,
        46, 159-176, `doi:10.1016/S0341-8162(01)00164-3`_.

    """

    pdalgo.pitparams(neighbours, area, diffnb, prop, propa, propb, fillTH, epsilon, boundPts)


def pit_stack(elev, allFill, sealevel):
    """
    This function calls a **pit filling algorithm** to compute depression-less elevation grid.

    Tip:
        Planchon and Darboux, 2001: A Fast, Simple and Versatile Algorithm to Fill the Depressions of Digital Elevation Models - Catena,
        46, 159-176, `doi:10.1016/S0341-8162(01)00164-3`_.

    .. _doi:10.1016/S0341-8162(01)00164-3: http://dx.doi.org/10.1016/S0341-8162(01)00164-3

    Args:
        elev: numpy arrays containing the nodes elevation.
        allFill: produce depression-less surface.
        sealevel: current elevation of sea level.

    Returns:
        - fillH - numpy array containing the filled elevations.


    Caution:
        The Planchon & Darboux (2001) algorithm is not as efficient as priority-queue approaches such as the one
        proposed in Barnes et al. (2014) and we now use this latest algorithm.

        Barnes, Lehman & Mulla 2014: An Efficient Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models -
        Computers & Geosciences, doi: 10.1016/`j.cageo.2013.01.009`_.

    .. _j.cageo.2013.01.009:  http://dx.doi.org/10.1016/j.cageo.2013.01.009

    """

    # Call stack based pit filling function from libUtils
    fillH = pdalgo.pitfilling(elev, allFill, sealevel)

    return fillH
