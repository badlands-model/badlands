subroutine SwanFindPoint ( x, y, kvert )
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
!   Finds the closest vertex index of the given point
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(out) :: kvert ! closest vertex index of given point
                                  ! Note: kvert = -1 indicate point is not found
    real, intent(in)     :: x     ! x-coordinate of given point
    real, intent(in)     :: y     ! y-coordinate of given point
!
!   Local variables
!
    integer, save                         :: ient = 0        ! number of entries in this subroutine
    integer                               :: iface           ! loop counter over faces
    integer                               :: ivert           ! loop counter over vertices
    integer                               :: v1              ! first vertex of present face
    integer                               :: v2              ! second vertex of present face
    !
    real                                  :: dismin          ! minimal distance found
    real                                  :: dist            ! computed distance
    real                                  :: dxb             ! x-component of length of boundary face
    real                                  :: dyb             ! y-component of length of boundary face
    real                                  :: r               ! relative distance of point to begin of boundary face
    real                                  :: reldis          ! relative distance of point to boundary face
    real                                  :: x1              ! x-coordinate of begin of boundary face
    real                                  :: x2              ! x-coordinate of end of boundary face
    real                                  :: xc              ! x-coordinate of closest vertex
    real                                  :: y1              ! y-coordinate of begin of boundary face
    real                                  :: y2              ! y-coordinate of end of boundary face
    real                                  :: yc              ! y-coordinate of closest vertex
    !
    logical                               :: SwanPointinMesh ! indicate whether a point is inside mesh
    !
    type(facetype), dimension(:), pointer :: face            ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanFindPoint')
    !
    ! point to face object
    !
    face => gridobject%face_grid
    !
    ! check whether point is outside the grid
    !
    if ( x < XCGMIN .or. x > XCGMAX .or. y < YCGMIN .or. y > YCGMAX ) then
       !
       kvert = -1
       return
       !
    endif
    !
    if ( SwanPointinMesh( x, y ) ) then
       !
       ! if point is inside the mesh then compute closest index
       !
       dismin = 1.e20
       !
       do ivert = 1, nverts
          !
          xc = xcugrd(ivert)
          yc = ycugrd(ivert)
          !
          dist = sqrt( (x-xc)**2 + (y-yc)**2 )
          if ( dist < dismin ) then
             kvert  = ivert
             dismin = dist
          endif
          !
       enddo
       !
    else
       !
       ! scan the boundary to look for the given point
       !
       ! loop over faces (both internal and boundary faces)
       !
       faceloop: do iface = 1, nfaces
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
             dxb = x2 - x1
             dyb = y2 - y1
             !
             reldis = abs( dyb*(x-x1) - dxb*(y-y1) ) / ( dxb*dxb + dyb*dyb )
             !
             if ( reldis < 0.01 ) then
                !
                r   = ( dxb*(x-x1) + dyb*(y-y1) ) / ( dxb*dxb + dyb*dyb )
                !
                if ( r < -0.01 .or. r > 1.01 ) then      ! 41.13
                   kvert = -1
                else
                   !
                   if ( r < 0.5 ) then
                      kvert = v1
                   else
                      kvert = v2
                   endif
                   exit faceloop
                   !
                endif
                !
             else
                kvert = -1
             endif
             !
          endif
          !
       enddo faceloop
       !
    endif
    !
end subroutine SwanFindPoint
