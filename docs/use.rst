###################
Usage
###################

.. image:: img/fig7.jpg
   :scale: 70 %
   :alt: strat
   :align: center

Portability
--------------

**Badlands** is an open-source package distributed under **GNU GPLv3 license**. The source code is available on `GitHub`_. Structurally the code is a *Python front end* with a *C* and *Fortran* middle layer to efficiently compute some of the heaviest functions. This Python-friendly version provides a programmable and flexible interface which maximises its portability across platforms.

.. _`GitHub`: https://github.com/badlands-model/badlands

.. tip::
  Instructions to install the code and the associated dependencies on a local system are provided in the next section of this documentation along with model options and a series of hands-on examples.

The easiest way to use **badlands** is via a *Docker container* (searching for *badlands-demo* on DockerHub or Kitematic) which is shipped with the complete list of dependencies, the model's companion, series of examples as well as workshop simulations.

.. caution::
  Models data and outputs ran from within the container will not persist when that container is no longer active.

To provide better interfacing between the container and the host filesystem, **badlands** image can be mounted on a local volume which allow for easy access and ability to store securely model results.

Interactions with other packages
--------------------------------

Triangulation
^^^^^^^^^^^^^

**Badlands** main calculations are performed on a TIN. However the code creates its own Delaunay triangulation using the Shewchuk’s Triangle library [1]_ from regularly defined input datasets (*e.g.*, topography grid, rain maps, tectonic maps). In that way, users provide standard regularly spaced *ASCII* datasets.

The only requirement is to follow a specific column-major order for the declaration of each nodes values which is consistent between imported datasets starting from the south-west and ending on the north-east corner.

Visualisation
^^^^^^^^^^^^^

Model results consist of time series of surface evolution, river and catchment dynamics grids as well as underlying stratigraphic architecture mesh. These outputs are all produced as **Hdf5 binary files** making it possible to interact with multiple existing visualisation and analysis software, such as Paraview_ or Visit_.

.. _Paraview: http://www.paraview.org
.. _Visit: https://wci.llnl.gov/simulation/computer-codes/visit/

Initial surface can be generated from UTM coordinates and functions have been added to easily extract Web Map Service dataset (one example is provided to illustrate how to define an initial topography grid from ETOPO5 datasets). Hdf5 files can also be quickly transformed in other conventional raster GIS file formats such as ASCII grids.

Flexural isostasy
^^^^^^^^^^^^^^^^^

To estimate flexural isostasy, *gFlex* modular python package [2]_ has been integrated as a component in **badlands**. It allows to compute isostatic deflections of Earth’s lithosphere with uniform or nonuniform flexural rigidity and couple the interactions with evolving surface loads induced by erosion and deposition associated to modelled surface processes.

.. Note::
      Flexural isostasy is obtained from the gflex_ package available on Github |:bomb:|

.. _gflex: https://github.com/awickert/gFlex

Geodynamic models
^^^^^^^^^^^^^^^^^

.. image:: img/collision_wedge.gif
   :scale: 100 %
   :alt: capability
   :align: center


uwgeodynamics_ [3]_ provides a way to couple an **underworld** model to **badlands**.

.. _uwgeodynamics: https://github.com/underworldcode/UWGeodynamics


.. code:: python

  import UWGeodynamics as GEO

  u = GEO.u
  air = GEO.Material()
  sediment = GEO.Material()
  Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(
    airIndex=[air.index], sedimentIndex=sediment.index,
    XML="resources/badlands.xml", resolution=1. * u.kilometre,
    checkpoint_interval=0.01 * u.megayears)

This will allow communication between the **uwgeodynamics** model and **badlands**. As in every **badlands** simulation, input parameters must be defined inside the *XML* file as described in this documentation.

.. role:: python(code)
    :language: py

The resulting :python:`Model` is a 2-way coupled thermo-mechanical model with surface processes, where the velocity field retrieved from the thermo-mechanical model is used to advect the surface in **badlands**. The surface in subject to both erosion and deposition. The distribution of materials in the thermo-mechanical (**underworld**) model is then updated.

.. important::
  It is recommended to use a higher spatial resolution in the surface processes model than in the thermo-mechanical model.


Users must define a list of material describing the air layers (usually, air and sticky air). It is also require to define an :python:`UWGeodynamics.Material` object describing the sediment that will be deposited. The index of the :python:`Material` is passed to the :python:`surfaceProcesses` function. Users can also provide an **underworld** function returning an index of an existing :python:`UWGeodynamics.Material`.


.. note::
  When the thermo-mechanical model is **2D**, the velocity field at the surface is extrapolated in the 3D dimension and the resulting model is a **T** or **2.5D** model (symmetric regional uplift). If the thermo-mechanical model is 3D the coupling is done in 3D.

----------

.. [1] J. R. Shewchuk -
  Triangle: Engineering a 2D quality mesh generator and Delaunay triangulator, pp. 203–222. Berlin, Heidelberg: Springer Berlin Heidelberg, 1996.

.. [2] A. D. Wickert -
  Open-source modular solutions for flexural isostasy: gflex v1.0, Geoscientific Model Development, vol. 9, no. 3, pp. 997–1017, 2016.


.. [3] Beucher et al. -
  UWGeodynamics: A teaching and research tool for numerical geodynamic modelling. Journal of Open Source Software, 4(36), 1136, `doi:10.21105/joss.01136`_, 2019.


.. _`doi:10.21105/joss.01136`: https://doi.org/10.21105/joss.01136
