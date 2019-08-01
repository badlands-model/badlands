subroutine SwanCheckGrid
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
!   Checks whether the grid is suited for computation
!
!   Method
!
!   For the following aspects the grid is checked for:
!       1) the number of cells around vertex must be at least 4 and not larger than 10
!       2) the angles in each triangle must be smaller than 143 degrees
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
!
    implicit none
!
!   Local variables
!
    integer                            :: i         ! loop counter
    integer, save                      :: ient = 0  ! number of entries in this subroutine
    integer                            :: j         ! loop counter
    integer                            :: v1        ! first vertex of present cell
    integer                            :: v2        ! second vertex of present cell
    integer                            :: v3        ! third vertex of present cell
    !
    real                               :: cosphi1   ! cosine of the angle of first vertex in a triangle
    real                               :: cosphi2   ! cosine of the angle of second vertex in a triangle
    real                               :: cosphi3   ! cosine of the angle of third vertex in a triangle
    real                               :: len12     ! squared length of face between vertices 1 and 2
    real                               :: len13     ! squared length of face between vertices 1 and 3
    real                               :: len23     ! squared length of face between vertices 2 and 3
    real                               :: xdif12    ! difference in x-coordinate of vertices 1 and 2
    real                               :: xdif13    ! difference in x-coordinate of vertices 1 and 3
    real                               :: xdif23    ! difference in x-coordinate of vertices 2 and 3
    real                               :: ydif12    ! difference in y-coordinate of vertices 1 and 2
    real                               :: ydif13    ! difference in y-coordinate of vertices 1 and 3
    real                               :: ydif23    ! difference in y-coordinate of vertices 2 and 3
    !
    integer, dimension(:), allocatable :: vcount    ! counts number of cells around each vertex
    !
    logical                            :: acute     ! indicates whether triangles are acute (.TRUE.) or not (.FALSE.)
    logical                            :: badvertex ! indicates vertex has too less cells surrounded
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanCheckGrid')
    !
    ! check whether the number of cells that meet at each internal vertex is at least 4
    ! (vertices at the boundaries not taken into account)
    ! check also whether the number of cells around each vertex is not larger than 10
    !
    allocate (vcount(nverts))
    vcount = 0
    !
    do i = 1, ncells
       !
       v1 = kvertc(1,i)
       v2 = kvertc(2,i)
       v3 = kvertc(3,i)
       !
       do j = 1, nverts
          if ( v1 == j ) vcount(j) = vcount(j) + 1
          if ( v2 == j ) vcount(j) = vcount(j) + 1
          if ( v3 == j ) vcount(j) = vcount(j) + 1
       enddo
       !
    enddo
    !
    badvertex = .false.
    do i = 1, nverts
       if ( vcount(i) > 0 .and. vcount(i) < 4 .and. vmark(i) == 0 ) badvertex = .true.
!ADC       if ( vcount(i) > 0 .and. vcount(i) < 4 .and. vmark(i) == 0 ) vmark(i) = 1
       if ( vcount(i) > 10 ) badvertex = .true.
    enddo
    !
    if ( badvertex ) call msgerr (2, 'number of cells around vertex is smaller than 4 or larger than 10')
    deallocate(vcount)
    !
    acute = .true.
    !
    ! check whether the angles in each triangle are smaller than a certain value (=143 degrees)
    !
    do i = 1, ncells
       !
       v1 = kvertc(1,i)
       v2 = kvertc(2,i)
       v3 = kvertc(3,i)
       !
       xdif12 = xcugrd(v2) - xcugrd(v1)
       ydif12 = ycugrd(v2) - ycugrd(v1)
       xdif13 = xcugrd(v3) - xcugrd(v1)
       ydif13 = ycugrd(v3) - ycugrd(v1)
       xdif23 = xcugrd(v3) - xcugrd(v2)
       ydif23 = ycugrd(v3) - ycugrd(v2)
       !
       len12 = xdif12*xdif12 + ydif12*ydif12
       len13 = xdif13*xdif13 + ydif13*ydif13
       len23 = xdif23*xdif23 + ydif23*ydif23
       !
       ! is triangle acute ?
       !
       if (acute) acute = (len12+len23>len13) .and. (len23+len13>len12) .and. (len13+len12>len23)
       !
       cosphi1 =( xdif12*xdif13 + ydif12*ydif13)/(sqrt(len12*len13))
       cosphi2 =(-xdif12*xdif23 - ydif12*ydif23)/(sqrt(len12*len23))
       cosphi3 =( xdif13*xdif23 + ydif13*ydif23)/(sqrt(len13*len23))
       !
       if ( cosphi1 <= -0.8 .or. cosphi2 <= -0.8 .or. cosphi3 <= -0.8 ) then
          call msgerr (2, 'an angle in a triangle is too large ')
       endif
       !
    enddo
    !
    if (acute) call msgerr (0, 'The grid contains solely acute triangles ')
    !
end subroutine SwanCheckGrid
