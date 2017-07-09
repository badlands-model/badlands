##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
   Top-level Badlands implementation.
"""

from .flow import flowNetwork
from .surface import FVmethod
from .surface import raster2TIN
from .simulation import oceanDyn
from .underland import eroMesh
from .underland import strataMesh
from .underland import stratiWedge
from .underland import carbMesh
from .flow import visualiseFlow
from .hillslope import diffLinear
from .surface import elevationTIN
from .surface import partitionTIN
from .surface import visualiseTIN
from .surface import visSurf
from .forcing import xmlParser
from .forcing import forceSim
from .forcing import isoFlex
from .forcing import carbGrowth
from .forcing import pelagicGrowth
from .simulation import buildMesh
from .simulation import checkPoints
from .simulation import buildFlux

# To build documentation, readthedocs.io loads pyBadlands into a virtualenv.
# It can't load the binary parts, though. To allow pyBadlands to be imported
# without the binaries, we check for the 'READTHEDOCS' environment variable.
import os
if 'READTHEDOCS' not in os.environ:
    from .libUtils import *
