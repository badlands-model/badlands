module SwanGridobjects
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
!   40.92, June 2008: new attribute for vertices: BPOL
!
!   Purpose
!
!   Module containing data structure for unstructured grids
!
!   Method
!
!   This module contains derived types for grid vertices, cells and faces
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
!
!   parameters for vertices
!
    integer, parameter :: MINVERTCELL = 4   ! mininum number of cells around a vertex
    integer, parameter :: MAXVERTCELL = 10  ! maximum number of cells around a vertex
    !
    integer, parameter :: MAXVERTATTI = 7   ! maximum number of attributes of type integer
                                            ! in data structure of vertices
    integer, parameter :: VERTID  = 1       ! identification number
    integer, parameter :: VMARKER = 2       ! boundary marker for vertices
                                            ! vert(j)%atti(VMARKER) = 1 is a boundary vertex
                                            ! vert(j)%atti(VMARKER) = 0 is an internal vertex
    integer, parameter :: VBC     = 3       ! marker for vertices with boundary condition
                                            ! vert(j)%atti(VBC) = 1 is a vertex where b.c. is given
                                            ! vert(j)%atti(VBC) = 0 is a vertex where no b.c. is given
    integer, parameter :: VERTF1  = 4       ! first face connecting to present vertex
    integer, parameter :: VERTF2  = 5       ! second face connecting to present vertex
    integer, parameter :: BINDX   = 6       ! indices for boundary points in ascending order
    integer, parameter :: BPOL    = 7       ! sequence number of boundary polygons
    !
    integer, parameter :: MAXVERTATTR = 2   ! maximum number of attributes of type real
                                            ! in data structure of vertices
    integer, parameter :: VERTX = 1         ! x-component of vertex
    integer, parameter :: VERTY = 2         ! y-component of vertex
!
!   parameters for cells
!
    integer, parameter :: MINCELLVERT = 3   ! mininum number of vertices in a cell
    integer, parameter :: MAXCELLVERT = 5   ! maximum number of vertices in a cell
    integer, parameter :: MINCELLFACE = 3   ! mininum number of faces in a cell
    integer, parameter :: MAXCELLFACE = 5   ! maximum number of faces in a cell
    !
    integer, parameter :: MAXCELLATTI = 7   ! maximum number of attributes of type integer
                                            ! in data structure of cells
    integer, parameter :: CELLID   = 1      ! identification number
    integer, parameter :: CMARKER  = 2      ! boundary marker for cells
                                            ! cell(j)%atti(CMARKER) = 1 is a boundary cell
                                            ! cell(j)%atti(CMARKER) = 0 is an internal cell
    integer, parameter :: NEXTCELL = 3      ! cell next to considered cell in counterclockwise direction
    integer, parameter :: CELLRECT = 4      ! /=0 if cell is rectangle
    integer, parameter :: CELLV1   = 5      ! first vertex number of present cell
    integer, parameter :: CELLV2   = 6      ! second vertex number of present cell
    integer, parameter :: CELLV3   = 7      ! third vertex number of present cell
    !
    integer, parameter :: MAXCELLATTR = 3   ! maximum number of attributes of type real
                                            ! in data structure of cells
    integer, parameter :: CELLAREA = 1      ! area of cell
    integer, parameter :: CELLCX   = 2      ! x-component of cell centroid
    integer, parameter :: CELLCY   = 3      ! y-component of cell centroid
    !
    ! parameters for faces
    !
    integer, parameter :: MAXFACEATTI = 6   ! maximum number of attributes of type integer
                                            ! in data structure of faces
    integer, parameter :: FACEID  = 1       ! identification number
    integer, parameter :: FMARKER = 2       ! boundary marker for faces
                                            ! face(j)%atti(FMARKER) = 1 is a boundary face
                                            ! face(j)%atti(FMARKER) = 0 is an internal face
    integer, parameter :: FACEV1  = 3       ! first vertex number of present face
    integer, parameter :: FACEV2  = 4       ! second vertex number of present face
    integer, parameter :: FACEC1  = 5       ! first cell number containing present face
    integer, parameter :: FACEC2  = 6       ! second cell number containing present face
                                            ! The ordening is such, that:
                                            !  - face(j)%atti(FACEC1) < face(j)%atti(FACEC2) always
                                            !  - the first cell lies left of translation vector
                                            !    (second vertex - first vertex) and the second cell
                                            !    consequently lies right of this vector
    !
    integer, parameter :: MAXFACEATTR = 6   ! maximum number of attributes of type real
                                            ! in data structure of faces
    integer, parameter :: FACELEN   = 1     ! length of face
    integer, parameter :: FACENORMX = 2     ! x-component of normal to present face
    integer, parameter :: FACENORMY = 3     ! y-component of normal to present face
    integer, parameter :: FACEMX    = 4     ! x-component of midpoint of present face
    integer, parameter :: FACEMY    = 5     ! y-component of midpoint of present face
    integer, parameter :: FACELINPF = 6     ! the interpolation factor when interpolating two quantities
                                            ! in cell centroids adjacent to face j according to
                                            ! q_face = q_m + face(j)%attr(FACELINPF) (q_n - q_m),
                                            ! where m = face(j)%atti(FACEC1) and n = face(j)%atti(FACEC2)
!
!   Module variables
!
!   ---
!
!   Source text
!
    type geomtype
      real :: det                               ! determinant
      real :: dx1 , dx2                         ! two components of covariant base vector a_(1)
      real :: dy1 , dy2                         ! two components of covariant base vector a_(2)
      real :: rdx1, rdx2                        ! first component of contravariant base vector rdx(b) = a^(b)_1
      real :: rdy1, rdy2                        ! second component of contravariant base vector rdy(b) = a^(b)_2
      real :: th1 , th2                         ! direction of faces pointing to present vertex
    end type geomtype

    type facetype
      integer :: atti(MAXFACEATTI)
      real    :: attr(MAXFACEATTR)
    end type facetype

    type celltype
      integer        :: nov                     ! actual number of vertices in a cell
      integer        :: nof                     ! actual number of faces in a cell
      logical        :: active                  ! true if cell is active
      type(geomtype) :: geom(MAXCELLVERT)
      type(facetype) :: face(MAXCELLFACE)
      integer        :: atti(MAXCELLATTI)
      real           :: attr(MAXCELLATTR)
      integer        :: updated(MAXCELLVERT)
    end type celltype

    type verttype
      integer        :: noc                     ! actual number of cells around a vertex
      logical        :: active                  ! true if vertex is active
      logical        :: fullupdated             ! true if vertex is updated in both geographic and spectral spaces
      type(celltype) :: cell(MAXVERTCELL)
      integer        :: atti(MAXVERTATTI)
      real           :: attr(MAXVERTATTR)
      integer        :: updated(MAXVERTCELL)
    end type verttype

    type gridtype
      type (verttype), dimension(:), pointer :: vert_grid
      type (celltype), dimension(:), pointer :: cell_grid
      type (facetype), dimension(:), pointer :: face_grid
    end type gridtype
    type(gridtype), target :: gridobject

end module SwanGridobjects
