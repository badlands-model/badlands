module SwanGriddata
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
!
!   Updates
!
!   40.80, July 2007: New Module
!
!   Purpose
!
!   Module containing data of unstructured grid
!
!   Method
!
!   Data with respect to unstructured grid need to be filled by a grid generator
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter :: meth_adcirc   = 1 !
    integer, parameter :: meth_triangle = 2 !
    integer, parameter :: meth_easy     = 3 !
!
!   Module variables
!
    integer                                    :: grid_generator ! used grid generator
    !
    integer                                    :: ncells         ! number of cells in (subdomain) grid
    integer                                    :: ncellsg        ! number of cells in global grid
    integer                                    :: nfaces         ! number of faces in grid
    integer                                    :: nverts         ! number of vertices in (subdomain) grid
    integer                                    :: nvertsg        ! number of vertices in global grid
    !
    integer, dimension(:,:), save, allocatable :: kvertc         !
    integer, dimension(:,:), save, allocatable :: kvertf         !
    !
    integer, dimension(:), save, allocatable   :: vmark          ! boundary marker for vertices
    !
    real                                       :: maxgsiz        ! maximum gridsize
    real                                       :: mingsiz        ! minimum gridsize
    !
    real, dimension(:), save, allocatable      :: xcugrd         ! the x-coordinates of the grid vertices
    real, dimension(:), save, allocatable      :: ycugrd         ! the y-coordinates of the grid vertices
!
!   Source text
!
end module SwanGriddata
