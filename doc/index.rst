pyBadlands
==========

Basin and Landscape Dynamics (Badlands) is a parallel TIN-based landscape evolution model, built to simulate topography development at various space and time scales. The model is capable of simulating hillslope processes (linear & non-linear diffusion), fluvial incision ('modified' Stream Power Law, Transport Capacity Law both for sediment erosion/transport/deposition), spatially and temporally varying geodynamic (horizontal + vertical displacements) and climatic forces which can be used to simulate changes in base level, as well as effects of climate changes or sea-level fluctuations. The model uses gFlex package which is designed to solve elastic plate flexure for applications to Earth's lithosphere.

Usage documentation and tutorials can be found on the .. _pyBadlands Wiki: https://github.com/badlands-model/pyBadlands/wiki

The current XML input file format is documented at .. _XML input for Badlands version 2: https://github.com/badlands-model/pyBadlands/wiki/XML-input-for-Badlands-version-2

Most users will only need to read the :doc:`public`

If you want to understand the internals of pyBadlands, you might be interested in the :doc:`internal`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   public
   internal

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
