pyBadlands
==========

Basin and Landscape Dynamics (Badlands) is a parallel TIN-based landscape evolution model, built to simulate topography development at various space and time scales. The model is capable of simulating hillslope processes (**linear** & **non-linear** diffusion), fluvial incision (*'modified'* **Stream Power Law**, **Transport Capacity Law** both for sediment  erosion/transport/deposition), spatially and temporally varying geodynamic (horizontal + vertical displacements) and climatic forces which can be used to simulate changes in base level, as well as effects of climate changes or sea-level fluctuations. The model uses `gFlex <https://github.com/awickert/gFlex>`_ package which is designed to solve elastic plate flexure for applications to Earth's lithosphere.

Usage documentation and tutorials can be found on the `pyBadlands Wiki <https://github.com/badlands-model/pyBadlands/wiki>`_

The current XML input file format is documented at `XML input for Badlands version 2 <https://github.com/badlands-model/pyBadlands/wiki/XML-input-for-Badlands-version-2>`_

Most users will only need to read the :doc:`models`

The remainder of the API documentation is interesting if you want to understand the internals of pyBadlands.

Contents
========

.. toctree::
   models
   flow
   forcing
   hillslope
   simulation
   surface
   underland

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
