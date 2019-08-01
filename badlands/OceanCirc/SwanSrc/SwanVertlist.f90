subroutine SwanVertlist
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
!   41.07: Casey Dietrich
!
!   Updates
!
!   40.80, July 2007: New subroutine
!   41.07, July 2009: small fix (assign ref.point to deepest point in case of no b.c.)
!
!   Purpose
!
!   Makes vertex list with specific order
!
!   Method
!
!   Sorting based on distance with respect to a reference point
!   This reference point can be either a vertex with boundary condition or a deepest point
!
!   Modules used
!
    use ocpcomm4
    use m_genarr, only: DEPTH
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Local variables
!
    integer, save                   :: ient = 0   ! number of entries in this subroutine
    integer                         :: istat      ! indicate status of allocation
    integer                         :: itmp       ! temporary stored integer for swapping
    integer                         :: j          ! loop counter over vertices
    integer                         :: k          ! counter
    integer, dimension(1)           :: kd         ! location of minimum value in array dist
    integer, dimension(1)           :: kx         ! location of minimum value in array of x-coordinates of boundary vertices
    integer, dimension(1)           :: ky         ! location of minimum value in array of y-coordinates of boundary vertices
    !
    real                            :: d1         ! distance of a point to origin
    real                            :: d2         ! distance of another point to origin
    real                            :: depmax     ! maximum depth found
    real                            :: rtmp       ! temporary stored real for swapping
    real                            :: x0         ! x-coordinate of reference point
    real                            :: y0         ! y-coordinate of reference point
    !
    real, dimension(:), allocatable :: dist       ! distance of each point with respect to reference point
    !
    type(verttype), dimension(:), pointer :: vert ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanVertlist')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    if(.not.allocated(vlist)) allocate (vlist(nverts), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanVertlist: array vlist ' )
       return
    endif
    !
    ! determine reference point
    !
    if ( all(mask=vert(:)%atti(VBC)==0) ) then
       !
       ! if no boundary condition is given then find the vertex with the maximum depth
       !
       depmax = -999.
       !
       do j = 1, nverts
          !
          if ( DEPTH(j) > depmax ) then
             !
             depmax = DEPTH(j)
             kx(1)  = j
             ky(1)  = j
             !
          endif
          !
       enddo
       !
    else
       !
       ! reference point is one of the vertices nearest to the origin where boundary condition is given
       !
       kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VBC)/=0)
       ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VBC)/=0)
       !
    endif
    !
    if ( kx(1) == ky(1) ) then
       x0 = vert(kx(1))%attr(VERTX)
       y0 = vert(ky(1))%attr(VERTY)
    else
       !
       d1 = sqrt((vert(kx(1))%attr(VERTX))**2+(vert(kx(1))%attr(VERTY))**2)
       d2 = sqrt((vert(ky(1))%attr(VERTX))**2+(vert(ky(1))%attr(VERTY))**2)
       !
       if ( d1 < d2 ) then
          x0 = vert(kx(1))%attr(VERTX)
          y0 = vert(kx(1))%attr(VERTY)
       else
          x0 = vert(ky(1))%attr(VERTX)
          y0 = vert(ky(1))%attr(VERTY)
       endif
       !
    endif
    !
    ! calculate distance of each point with respect to reference point
    !
    allocate (dist(nverts))
    do j = 1, nverts
       dist(j) = sqrt((vert(j)%attr(VERTX)-x0)**2 + (vert(j)%attr(VERTY)-y0)**2)
    enddo
    !
    ! sort vertex list in order of increasing distance
    !
    vlist=(/ (j, j=1, nverts) /)
    !
    do j = 1, nverts-1
       !
       kd = minloc(dist(j:nverts))
       k  = kd(1) + j-1
       !
       if ( k /= j ) then
          !
          rtmp     = dist(j)
          dist(j)  = dist(k)
          dist(k)  = rtmp
          !
          itmp     = vlist(j)
          vlist(j) = vlist(k)
          vlist(k) = itmp
          !
       endif
       !
    enddo
    !
    deallocate(dist)
    !
end subroutine SwanVertlist
