subroutine SwanGridVert ( nverts, xcugrd, ycugrd, vmark )
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
!   Fills vertex-based data structure
!
!   Method
!
!   Based on unstructured grid
!
!   Modules used
!
    use ocpcomm4
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                       :: nverts  ! number of vertices in grid
    !
    integer, dimension(nverts), intent(in)    :: vmark   ! boundary marker for vertices
    !
    real, dimension(nverts), intent(in)       :: xcugrd  ! the x-coordinates of the grid vertices
    real, dimension(nverts), intent(in)       :: ycugrd  ! the y-coordinates of the grid vertices
!
!   Local variables
!
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: ivert    ! loop counter over vertices
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGridVert')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    ! loop over all vertices
    !
    do ivert = 1, nverts
       !
       ! identification number
       !
       vert(ivert)%atti(VERTID) = ivert
       !
       ! marks boundary vertex
       !
       vert(ivert)%atti(VMARKER) = min(1,vmark(ivert))
       !
       ! initially, this vertex is no boundary condition point
       !
       vert(ivert)%atti(VBC) = 0
       !
       ! initially, this vertex is active
       !
       vert(ivert)%active = .true.
       !
       ! store coordinates
       !
       vert(ivert)%attr(VERTX) = xcugrd(ivert)
       vert(ivert)%attr(VERTY) = ycugrd(ivert)
       !
    enddo

end subroutine SwanGridVert
