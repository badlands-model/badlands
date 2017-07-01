module SwanCompdata
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2011  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.80: Marcel Zijlema
!   40.92: Marcel Zijlema
!
!   Updates
!
!   40.80, July 2007: New Module
!   40.92, June 2008: changes with respect to boundary polygons
!
!   Purpose
!
!   Module containing data for computation with unstructured grid
!
!   Method
!
!   Data based on unstructured grid
!
!   Modules used
!
    use swcomm3
!
    implicit none
!
!   Module parameters
!
!
!   Module variables
!
    integer                                    :: nbpol  ! total number of boundary polygons
    integer, dimension(10000)                  :: nbpt   ! number of boundary vertices for each boundary polygon
!
    integer, dimension(MICMAX)                 :: vs     ! computational stencil, i.e. set of vertices
                                                         ! needed for the computation of a new value
                                                         ! in the present vertex
!$omp threadprivate(vs)
!
    integer, dimension(:,:), save, allocatable :: blist  ! list of boundary vertices in ascending order for each boundary polygon
    integer, dimension(:)  , save, allocatable :: vlist  ! vertex list
!
!   Source text
!
end module SwanCompdata
