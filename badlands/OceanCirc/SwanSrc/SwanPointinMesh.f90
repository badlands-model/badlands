logical function SwanPointinMesh ( x, y )
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
!   40.80, June 2007: New subroutine
!
!   Purpose
!
!   Checks whether the given point is inside the mesh
!
!   Method
!
!   draw a vertical line from the point and count the number of crossings with boundary faces
!   if the number of crossings is odd then the given point is inside the mesh
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, intent(in) :: x ! x-coordinate of given point
    real, intent(in) :: y ! y-coordinate of given point
!
!   Local variables
!
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! loop counter over faces
    integer                               :: numcrs   ! number of crossings with boundary faces
    integer                               :: v1       ! first vertex of present face
    integer                               :: v2       ! second vertex of present face
    !
    real                                  :: x1       ! x-coordinate of begin of boundary face
    real                                  :: x2       ! x-coordinate of end of boundary face
    real                                  :: y1       ! y-coordinate of begin of boundary face
    real                                  :: y2       ! y-coordinate of end of boundary face
    real                                  :: yc       ! y-coordinate of cross point
    !
    type(facetype), dimension(:), pointer :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanPointinMesh')
    !
    ! point to face object
    !
    face => gridobject%face_grid
    !
    numcrs = 0
    !
    ! loop over faces (both internal and boundary faces)
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 1 ) then   ! boundary face
          !
          v1 = face(iface)%atti(FACEV1)
          v2 = face(iface)%atti(FACEV2)
          !
          x1 = xcugrd(v1)
          y1 = ycugrd(v1)
          x2 = xcugrd(v2)
          y2 = ycugrd(v2)
          !
          if ( ( (x1 > x) .and. (x2 <= x) ) .or. ( (x2 > x) .and. (x1 <= x) ) ) then
             !
             if ( y1 > y .or. y2 > y ) then
                !
                yc = y1 + (x-x1) * (y2-y1) / (x2-x1)
                if ( yc > y ) numcrs = numcrs + 1
                !
             endif
             !
          endif
          !
       endif
       !
    enddo
    !
    ! if number of crossings is odd then point is inside the grid
    !
    if ( mod(numcrs,2) == 1 ) then
       SwanPointinMesh = .true.
    else
       SwanPointinMesh = .false.
    endif
    !
end function SwanPointinMesh
