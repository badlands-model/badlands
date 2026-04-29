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
import meshplex
import warnings
import triangle as triangle

import os

if "READTHEDOCS" not in os.environ:
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
        self.maxNgbh = None
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
        walltime = time.process_time()

        node_coords = numpy.zeros((len(self.node_coords), 3))
        node_coords[:, :2] = self.node_coords
        Tmesh = meshplex.MeshTri(node_coords, self.cells[:, :3])
        self.control_volumes = numpy.abs(Tmesh.control_volumes)
        self.control_volumes[numpy.isnan(self.control_volumes)] = 1.0

        # Voronoi and simplices declaration
        Tmesh.create_edges()
        cc = Tmesh.cell_circumcenters
        if meshplex.__version__ >= "0.16.0":
            edges_nodes = Tmesh.edges["points"]
            cells_nodes = Tmesh.cells("points")
            cells_edges = Tmesh.cells("edges")
        elif meshplex.__version__ >= "0.14.0":
            edges_nodes = Tmesh.edges["points"]
            cells_nodes = Tmesh.cells["points"]
            cells_edges = Tmesh.cells["edges"]
        else:
            edges_nodes = Tmesh.edges["nodes"]
            cells_nodes = Tmesh.cells["nodes"]
            cells_edges = Tmesh.cells["edges"]

        # Finite volume discretisation
        (
            self.neighbours,
            self.vor_edges,
            self.edge_length,
            maxNgbhs,
        ) = fvframe.definetin(node_coords, cells_nodes, cells_edges, edges_nodes, cc.T,)
        if verbose:
            print(
                " - construct Finite Volume representation ",
                time.process_time() - walltime,
            )

        # Maximum number of neighbours for each partition
        self.maxNgbh = numpy.array(maxNgbhs)

        return

    def construct_FV(self, lGIDs, verbose=False):
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
            lGIDs: numpy integer-type array filled with the global vertex IDs for each local grid located within the partition (including those on the edges).
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).


        Important:
            In the background the code is calling **Triangle** a python wrapper around *Jonathan Richard
            Shewchuk*'s two-dimensional quality mesh generator and delaunay triangulator library,
            available `here <quake_>`_.  The source is available on Github_.

        .. _quake: http://www.cs.cmu.edu/~quake/triangle.html
        .. _Github: https://github.com/drufat/triangle
        """

        # Call finite volume function
        self._FV_utils(lGIDs)

        return
