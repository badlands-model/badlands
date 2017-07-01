subroutine SwanFindObstacles ( cross )
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
!   40.80, March 2008: New subroutine
!
!   Purpose
!
!   Search for obstacles in computational grid and store them
!
!   Method
!
!   for each face an obstacle is found when they crossed each other
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use m_obsta
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, dimension(nfaces), intent(out) :: cross           ! contains sequence number of obstacles for each face
                                                               ! where they crossing or zero if no crossing
!
!   Local variables
!
    integer, save                         :: ient = 0          ! number of entries in this subroutine
    integer                               :: iface             ! loop counter over faces
    integer                               :: j                 ! loop counter
    integer                               :: k                 ! loop counter
    integer                               :: vb                ! vertex of begin of present face
    integer                               :: ve                ! vertex of end of present face
    !
    real, dimension(2)                    :: xobs              ! x-coordinate of obstacle point
    real, dimension(2)                    :: xv                ! x-coordinate of vertex of face
    real, dimension(2)                    :: yobs              ! y-coordinate of obstacle point
    real, dimension(2)                    :: yv                ! y-coordinate of vertex of face
    !
    logical                               :: SwanCrossObstacle ! indicate whether a face cross an obstacle
    !
    type(OBSTDAT), pointer                :: cobst             ! pointer to a considered obstacle
    !
    type(facetype), dimension(:), pointer :: face              ! datastructure for faces with their attributes
    type(verttype), dimension(:), pointer :: vert              ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanFindObstacles')
    !
    ! point to vertex and face objects
    !
    vert => gridobject%vert_grid
    face => gridobject%face_grid
    !
    ! initialize array cross
    !
    cross = 0
    !
    ! go to first obstacle
    !
    cobst => FOBSTAC
    !
    do j = 1, NUMOBS
       !
       xobs(1) = cobst%XCRP(1)
       yobs(1) = cobst%YCRP(1)
       !
       do k = 2, cobst%NCRPTS    ! number of corner points of considered obstacle
          !
          xobs(2) = cobst%XCRP(k)
          yobs(2) = cobst%YCRP(k)
          !
          ! loop over faces
          !
          faceloop : do iface = 1, nfaces
             !
             if ( face(iface)%atti(FMARKER) == 1 ) cycle faceloop    ! boundary face
             !
             vb = face(iface)%atti(FACEV1)
             ve = face(iface)%atti(FACEV2)
             !
             xv(1) = vert(vb)%attr(VERTX)
             yv(1) = vert(vb)%attr(VERTY)
             xv(2) = vert(ve)%attr(VERTX)
             yv(2) = vert(ve)%attr(VERTY)
             !
             if ( SwanCrossObstacle( xv, yv, xobs, yobs ) ) cross(iface) = j
             !
          enddo faceloop
          !
          xobs(1) = xobs(2)
          yobs(1) = yobs(2)
          !
       enddo
       !
       ! go to next obstacle, if present
       !
       if (.not.associated(cobst%NEXTOBST)) exit
       cobst => cobst%NEXTOBST
       !
    enddo
    !
end subroutine SwanFindObstacles
