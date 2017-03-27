pyBadlands - python version of Badlands
=====
    
<div align="center">
    <img width=1000 src="https://github.com/badlands-model/Badlands-doc/blob/master/figures/StratigraphciModelling.png" alt="sketch Badlands" title="sketch of Badlands range of models."</img>
</div>

[![DOI](https://zenodo.org/badge/51286954.svg)](https://zenodo.org/badge/latestdoi/51286954)

## Overview

**Ba**sin an**d** **Lan**dscape **D**ynamic**s** (**Badlands**) is a parallel TIN-based landscape evolution model, built to simulate topography development at various space and time scales. The model is capable of simulating hillslope processes (**linear** & **non-linear** diffusion), fluvial incision (*'modified'* **Stream Power Law**, **Transport Capacity Law** both for sediment  erosion/transport/deposition), spatially and temporally varying geodynamic (horizontal + vertical displacements) and climatic forces which can be used to simulate changes in base level, as well as effects of climate changes or sea-level fluctuations. The model uses  [**gFlex**](https://github.com/awickert/gFlex) package which is designed to solve elastic plate flexure for applications to Earth's lithosphere.

[Badlands overview](https://prezi.com/5y1usorz8e8k/badlands-overview/?utm_campaign=share&utm_medium=copy) - Basin Genesis Hub presentation (2017)

## Getting started

For installation information and documentation visit our github [**wiki page**](https://github.com/badlands-model/pyBadlands/wiki) which provides several useful notes on how to start using the tool.

The easiest way to get started is with the [Docker container](https://hub.docker.com/u/badlandsmodel/) using [Kitematic](https://docs.docker.com/kitematic/userguide/). Once **Kitematic** is installed on your computer, open it and look for **pybadlands-demo** via the *search* menu.

If you want to install it yourself, these 2 Dockerfiles ([**dependencies**](https://github.com/badlands-model/pyBadlands-Dependencies-Docker/blob/master/Dockerfile) & [**code**](https://github.com/badlands-model/pyBadlands-Docker-Demo/blob/master/Dockerfile)) are the best documentation of the required packages.

The latest pyBadlands version is the one that’s in our Github [repository](https://github.com/badlands-model/pyBadlands). Get it using this shell command, which requires Git: 
* `git clone https://github.com/badlands-model/pyBadlands.git`

**API documentation** is available from [http://pybadlands.readthedocs.io](http://pybadlands.readthedocs.io/en/master/)

## The specs...

The model is based on the following characteristics:
* The finite volume approach from [**Tucker et al. (2001)**](http://www.sciencedirect.com/science/article/pii/S0098300400001345) based on the dual Delaunay-Voronoi framework is used to solve the continuity equation explicitly,
* Node ordering is perform efficiently based on the work from [**Braun & Willett (2013)**](http://www.sciencedirect.com/science/article/pii/S0169555X12004618),
* 3D surface deformations using the node refinement technique proposed by [**Thieulot et al. ( 2014)**](http://onlinelibrary.wiley.com/doi/10.1002/2014GC005490/abstract;jsessionid=48A885F79A40B1E3E76AFC1BEAA2B238.f03t03).
* Orographic precipitation using [**Smith & Barstad (2004)**](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(2004)061%3C1377%3AALTOOP%3E2.0.CO%3B2) linear model to compute topographic induced rain field.
* Varying erodibility layers (both horizontally and vertically) to simulate impact of changing sediment characteristics on landscape evolution. 

A set of functions for _pre_ & _post_-processing of **Badlands** inputs and outputs is available in a GitHub [**Companion**](https://github.com/badlands-model/pyBadlands-Companion) repository which is already shipped with the Badlands Docker container.

### Community driven

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.

### Versioning  

Current releases (**DOI** citable): 

_Release v1.0.0 :_  
+ **[Release 1: badlands-model/pyBadlands](http://doi.org/10.5281/zenodo.160412), Zenodo, doi:10.5281/zenodo.160412, 2016.**

_Release v1.1.0 :_  
+ **[Release 2: badlands-model/pyBadlands](http://doi.org/10.5281/zenodo.268053), Zenodo, doi:10.5281/zenodo.268053, 2017.**

_v2.0.0 Release Candidate :_
+ The ```master``` branch is a stable version of the code with the most recent developments embedded and is expected to be released as **version 2.0.0**

### Reporting  

If you come accross a bug or if you need some help compiling or using the code you can: 

- go through our mailing list [Archive](http://mailman.sydney.edu.au/pipermail/badlands/)
- subscribe to our [mailing list](http://mailman.sydney.edu.au/mailman/listinfo/badlands)
- or [drop us a line](mailto:badlands@mailman.sydney.edu.au)

## Hands-on examples

A compilation of notebooks with examples are proposed to give you a quick feeling of what could be done with the code. Testing models:

+ delta evolution under sea-level fluctuation [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/delta/delta.ipynb)
+ impact of climate on mountain dynamic [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/mountain/mountain.ipynb)
+ basin filling associated to a strike-slip fault system [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/strikeslip/strike-slip.ipynb)
+ infilling of a crater-type topography [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/crater/crater.ipynb)
+ flexural response due to loading and unloading under variable elastic thickness [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/flexure/flexure.ipynb)
+ quick setup of real topography/bathymetry model using etopo1 [nbviewer](http://nbviewer.jupyter.org/github/badlands-model/pyBadlands/blob/master/Examples/etopo/etopo.ipynb)

### Documentation related to the physics & assumptions

+ **Salles, T. & Hardiman, L.: [Badlands: An open-source, flexible and parallel framework to study landscape dynamics](http://dx.doi.org/10.1016/j.cageo.2016.03.011), Computers & Geosciences, 91, 77-89, doi:10.1016/j.cageo.2016.03.011, 2016.**
+ **Salles, T.: [Badlands: A parallel basin and landscape dynamics model](http://dx.doi.org/10.1016/j.softx.2016.08.005), SoftwareX, 5, 195–202, doi:10.1016/j.softx.2016.08.005, 2016.**

### Published research studies using Badlands:

+ **Salles, T., N. Flament, and D. Muller: [Influence of mantle flow on the drainage of eastern Australia since the Jurassic Period](http://dx.doi.org/10.1002/2016GC006617), Geochem. Geophys. Geosyst., 18, doi:10.1002/2016GC006617, 2017** -- [Supplementary materials: Australian Landscape Dynamic](https://github.com/badlands-model/g-cubed-2016)


When you use **Badlands**, please cite the above papers.
