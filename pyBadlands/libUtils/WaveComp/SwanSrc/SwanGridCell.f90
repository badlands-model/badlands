subroutine SwanGridCell ( ncells, nverts, xcugrd, ycugrd, kvertc )
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
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   Fills cell-based data structure
!
!   Method
!
!   Based on unstructured grid
!   Note: we restrict ourselves to triangles only!
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                       :: ncells  ! number of cells in grid
    integer, intent(in)                       :: nverts  ! number of vertices in grid
    !
    integer, dimension(3, ncells), intent(in) :: kvertc  ! vertices of the cell
                                                         ! (must be filled by a gridgenerator!)
    !
    real, dimension(nverts), intent(in)       :: xcugrd  ! the x-coordinates of the grid vertices
    real, dimension(nverts), intent(in)       :: ycugrd  ! the y-coordinates of the grid vertices
!
!   Local variables
!
    integer                               :: icell    ! loop counter over cells / index of present cell
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: ivert    ! loop counter over vertices / index of present vertex
    integer                               :: j        ! loop counter
    integer                               :: jcell    ! index of next cell
    integer                               :: k        ! auxiliary integer / loop counter
    integer                               :: l        ! loop counter
    integer                               :: noc      ! number of cells around considered vertex
    integer, parameter                    :: nov = 3  ! number of vertices in present cell (triangles only)
    integer, dimension(3)                 :: v        ! vertices in present cell
    integer                               :: vc       ! considered vertex
    integer                               :: vn       ! first upwave vertex of next cell
    integer                               :: vp       ! last upwave vertex of present cell
    !
    integer, dimension(:), allocatable    :: ivlist   ! list of index vertices
    !
    real                                  :: carea    ! area of the present cell (triangles only)
    real                                  :: dx1      ! first component of covariant base vector a_(1)
    real                                  :: dx2      ! second component of covariant base vector a_(1)
    real                                  :: dy1      ! first component of covariant base vector a_(2)
    real                                  :: dy2      ! second component of covariant base vector a_(2)
    real                                  :: rdet     ! reciproke of determinant
    real                                  :: th1      ! direction of one face pointing to present vertex
    real                                  :: th2      ! direction of another face pointing to present vertex
    real                                  :: thdiff   ! difference between th2 and th1
    real, dimension(2)                    :: vec12    ! translation vector of coordinates: vertex2 - vertex1
    real, dimension(2)                    :: vec13    ! translation vector of coordinates: vertex3 - vertex1
    real                                  :: xc       ! x-coordinate of the cell-centroid
    real                                  :: yc       ! y-coordinate of the cell-centroid
    !
    logical                               :: nxtcell  ! indicate whether there is next cell in counterclockwise direction
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGridCell')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! loop over all cells
    !
    do icell = 1, ncells
       !
       ! determine number of vertices and faces (triangles only!)
       !
       cell(icell)%nov = nov
       cell(icell)%nof = nov
       !
       ! identification number
       !
       cell(icell)%atti(CELLID) = icell
       !
       ! cell is triangle and initiatively active
       !
       cell(icell)%atti(CELLRECT) = 0
       cell(icell)%active         = .true.
       !
       ! store vertices of the cell
       !
       v(1) = kvertc(1,icell)
       v(2) = kvertc(2,icell)
       v(3) = kvertc(3,icell)
       cell(icell)%atti(CELLV1) = v(1)
       cell(icell)%atti(CELLV2) = v(2)
       cell(icell)%atti(CELLV3) = v(3)
       !
       ! store area of the cell in carea
       !
       vec12(1) = xcugrd(v(2)) - xcugrd(v(1))
       vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
       vec13(1) = xcugrd(v(3)) - xcugrd(v(1))
       vec13(2) = ycugrd(v(3)) - ycugrd(v(1))
       carea = 0.5*abs(vec12(1)*vec13(2) - vec13(1)*vec12(2))
       !
       ! store local covariant and contravariant base vectors at each vertex
       ! next, store directions of faces pointing to present vertex
       !
       dx1 = -vec12(1)
       dy1 = -vec12(2)
       dx2 = -vec13(1)
       dy2 = -vec13(2)
       !
       do j = 1, nov
          !
          rdet = 1./(dy2*dx1 - dy1*dx2)
          !
          cell(icell)%geom(j)%det  =  dy2*dx1 - dy1*dx2
          cell(icell)%geom(j)%dx1  =  dx1
          cell(icell)%geom(j)%dy1  =  dy1
          cell(icell)%geom(j)%dx2  =  dx2
          cell(icell)%geom(j)%dy2  =  dy2
          cell(icell)%geom(j)%rdx1 =  dy2*rdet
          cell(icell)%geom(j)%rdy1 = -dx2*rdet
          cell(icell)%geom(j)%rdx2 = -dy1*rdet
          cell(icell)%geom(j)%rdy2 =  dx1*rdet
          !
          th1 = atan2(dy1,dx1)
          th2 = atan2(dy2,dx2)
          !
          thdiff = th1 - th2
          do
             if ( abs(thdiff) <= PI ) exit
             th1 = th1 - sign (2., thdiff) * PI
             thdiff = th1 - th2
          enddo
          !
          cell(icell)%geom(j)%th1 = th1
          cell(icell)%geom(j)%th2 = th2
          !
          dx1 = dx2 - dx1
          dy1 = dy2 - dy1
          dx2 = dx1 - dx2
          dy2 = dy1 - dy2
          !
       enddo
       !
       ! determine orientation of the mesh
       !
       if ( icell == 1 ) then
          !
          if ( dy2*dx1 > dy1*dx2 ) then
             CVLEFT = .false.             ! right-handed
          else
             CVLEFT = .true.              ! left-handed
          endif
          !
       endif
       !
       ! store coordinates of centroid
       !
       xc = 0.
       yc = 0.
       do j = 1, nov
          xc = xc + xcugrd(kvertc(j,icell))
          yc = yc + ycugrd(kvertc(j,icell))
       enddo
       xc = xc / real(nov)
       yc = yc / real(nov)
       !
       cell(icell)%attr(CELLAREA) = carea
       cell(icell)%attr(CELLCX  ) = xc
       cell(icell)%attr(CELLCY  ) = yc
       !
       ! store vertices of each face of the cell
       !
       do j = 1, nov
          k = mod(j,nov)+1
          cell(icell)%face(j)%atti(FACEV1) = kvertc(j,icell)
          cell(icell)%face(j)%atti(FACEV2) = kvertc(k,icell)
       enddo
       !
    enddo
    !
    allocate(ivlist(nverts))
    !
    ! loop over all vertices
    !
    do ivert = 1, nverts
       !
       ! identify the considered vertex and store index
       !
       vc = vert(ivert)%atti(VERTID)
       ivlist(vc) = ivert
       !
       ! initialize number of cells around vertex
       !
       vert(ivert)%noc = 0
       !
    enddo
    !
    do icell = 1, ncells
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! add present cell to each of these vertices
       !
       do j = 1, nov
          !
          ivert = ivlist(v(j))
          noc   = vert(ivert)%noc +1
          !
          vert(ivert)%noc = noc
          vert(ivert)%cell(noc)%atti(CELLID) = icell
          !
       enddo
       !
    enddo
    !
    deallocate(ivlist)
    !
    do ivert = 1, nverts
       !
       noc = vert(ivert)%noc
       !
       ! loop over cells around considered vertex
       !
       do j = 1, noc
          !
          ! get cell and its vertices
          !
          icell = vert(ivert)%cell(j)%atti(CELLID)
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          ! pick up last upwave vertex (counterclockwise counting of vertices is assumed)
          !
          do k = 1, nov
             if ( v(k) == ivert ) then
                vp = v(mod(k+1,nov)+1)
                exit
             endif
          enddo
          !
          ! search for next cell in counterclockwise direction
          !
          nxtcell = .false.
          !
          do l = 1, noc
             !
             ! get a cell and its vertices
             !
             jcell = vert(ivert)%cell(l)%atti(CELLID)
             !
             v(1) = cell(jcell)%atti(CELLV1)
             v(2) = cell(jcell)%atti(CELLV2)
             v(3) = cell(jcell)%atti(CELLV3)
             !
             ! pick up first upwave vertex (counterclockwise counting of vertices is assumed)
             !
             do k = 1, nov
                if ( v(k) == ivert ) then
                   vn = v(mod(k,nov)+1)
                   exit
                endif
             enddo
             !
             ! check whether first upwave vertex of next cell equals last upwave vertex of present cell
             !
             if ( vn == vp ) then
                nxtcell = .true.
                exit
             endif
             !
          enddo
          !
          if ( .not.nxtcell ) jcell = 0
          vert(ivert)%cell(j)%atti(NEXTCELL) = jcell
          !
       enddo
       !
    enddo

end subroutine SwanGridCell
