---
layout: post
title: Introducing Badlands
---

![Example of Landscape evolution with Badlands](/images/riseofthephoenix.png)

## Overview

**Ba**sin an**d** **Lan**dscape **D**ynamic**s** (**Badlands**) is a parallel TIN-based landscape evolution model, built to simulate topography development at various space and time scales. The model is presently capable of simulating hillslope processes (**linear** diffusion), fluvial incision (*'modified'* **SPL**:  erosion/transport/deposition), spatially and temporally varying geodynamic (horizontal + vertical displacements) and climatic forces which can be used to simulate changes in base level, as well as effects of climate changes or sea-level fluctuations.

## The specs...

The model is mainly written in fortran and is based on the following characteristics:

* The finite volume approach from Tucker et al. (2001) based on the dual Delaunay-Voronoi framework is used to solve the continuity equation explicitly,
* Node ordering is perform efficiently based on the work from Braun & Willett (2013),
* A Hilbert Space-Filling Curve method algorithm (Zoltan) is used to partition the TIN-based surface into subdomains,
* Drainage network partitioning is generated through METIS library.

## Community driven

The code is conceived as an open-source project, and is an ideal tool for both **Research** and **Learning** purposes.
