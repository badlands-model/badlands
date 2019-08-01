##################
API Documentation
##################

**Badlands** is a long-term surface evolution model built to simulate landscape development, sediment transport and
sedimentary basins formation from upstream regions down to marine environments.

.. image:: img/component.png
   :scale: 30 %
   :alt: capability
   :align: center

The above figure illustrates the relationships between the different modules of the framework.

Below you will find the *nitty-gritty components* of the Python code...
The following API documentation provides a nearly complete description of the different functions
that compose the code... |:v:|

.. warning::
  It is worth mentioning that a set of **fortran** & **C** functions are also part of **badlands** code but are not described
  in this API...


Model class
-----------

.. automodule:: model
    :members:

Flow Network
------------

.. automodule:: flow
    :members:

flowNetwork
^^^^^^^^^^^

.. automodule:: flow.flowNetwork
    :members:

visualiseFlow
^^^^^^^^^^^^^

.. automodule:: flow.visualiseFlow
    :members:

Forcing
-------

.. automodule:: forcing
    :members:

forceSim
^^^^^^^^^^^^^

.. automodule:: forcing.forceSim
    :members:

carbGrowth
^^^^^^^^^^^^^

.. automodule:: forcing.carbGrowth
    :members:

isoFlex
^^^^^^^^^^^^^

.. automodule:: forcing.isoFlex
    :members:

pelagicGrowth
^^^^^^^^^^^^^

.. automodule:: forcing.pelagicGrowth
    :members:

xmlParser
^^^^^^^^^^^^^

.. automodule:: forcing.xmlParser
    :members:

Hillslope
---------

.. automodule:: hillslope
    :members:

diffLinear
^^^^^^^^^^

.. automodule:: hillslope.diffLinear
    :members:

Simulation
------------

.. automodule:: simulation
    :members:

buildFlux
^^^^^^^^^^

.. automodule:: simulation.buildFlux
    :members:

buildMesh
^^^^^^^^^^

.. automodule:: simulation.buildMesh
    :members:

checkPoints
^^^^^^^^^^^^

.. automodule:: simulation.checkPoints
    :members:

waveSed
^^^^^^^^^^

.. automodule:: simulation.waveSed
    :members:

Surface
--------

.. automodule:: surface
    :members:

FVmethod
^^^^^^^^^^^^

.. automodule:: surface.FVmethod
    :members:

elevationTIN
^^^^^^^^^^^^

.. automodule:: surface.elevationTIN
    :members:

partitionTIN
^^^^^^^^^^^^

.. automodule:: surface.partitionTIN
    :members:

raster2TIN
^^^^^^^^^^^^

.. automodule:: surface.raster2TIN
    :members:

visualiseTIN
^^^^^^^^^^^^

.. automodule:: surface.visualiseTIN
    :members:

Underland
---------

.. automodule:: underland
    :members:

carbMesh
^^^^^^^^^^^^

.. automodule:: underland.carbMesh
    :members:

eroMesh
^^^^^^^^^^^^

.. automodule:: underland.eroMesh
    :members:

strataMesh
^^^^^^^^^^^^

.. automodule:: underland.strataMesh
    :members:

stratiWedge
^^^^^^^^^^^^

.. automodule:: underland.stratiWedge
    :members:
