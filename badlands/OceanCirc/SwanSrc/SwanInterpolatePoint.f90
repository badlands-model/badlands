subroutine SwanInterpolatePoint ( foutp, x, y, finp, excval )
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
!   40.80, August 2007: New subroutine
!
!   Purpose
!
!   Interpolates given scalar to given point
!
!   Method
!
!   First, look for closest vertex and next, interpolate given scalar inside triangle where given point is resided
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, intent(in)                    :: excval ! exception value for given scalar
    real, dimension(nverts), intent(in) :: finp   ! given scalars defined on the computational grid
    real, intent(out)                   :: foutp  ! output scalar at given point
    real, intent(in)                    :: x      ! x-coordinate of given point
    real, intent(in)                    :: y      ! y-coordinate of given point
!
!   Local variables
!
    integer                               :: icell     ! cell index
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: ivert     ! vertex index
    integer                               :: jc        ! loop counter
    integer                               :: k         ! loop counter
    integer, dimension(3)                 :: v         ! vertices in present cell
    !
    real                                  :: carea     ! area of the present cell
    real                                  :: dxp       ! distance between given point and present vertex in x-direction
    real                                  :: dyp       ! distance between given point and present vertex in y-direction
    real                                  :: eps       ! a small number
    real                                  :: phi1      ! value of given scalar in first vertex of considered cell
    real                                  :: phi2      ! value of given scalar in second vertex of considered cell
    real                                  :: phi3      ! value of given scalar in third vertex of considered cell
    real                                  :: phic      ! value of given scalar in centroid of considered cell
    real                                  :: th        ! direction of given point to present vertex
    real                                  :: th1       ! direction of one face pointing to present vertex
    real                                  :: th2       ! direction of another face pointing to present vertex
    real                                  :: thdiff    ! difference between th and th2
    real, dimension(2)                    :: vec12     ! translation vector of coordinates: vertex2 - vertex1
    real, dimension(2)                    :: vec23     ! translation vector of coordinates: vertex3 - vertex2
    real, dimension(2)                    :: vec31     ! translation vector of coordinates: vertex1 - vertex3
    real                                  :: xc        ! x-coordinate of the cell-centroid
    real                                  :: yc        ! y-coordinate of the cell-centroid
    real                                  :: xgrs      ! x-component of gradient scalar vector
    real                                  :: ygrs      ! y-component of gradient scalar vector
    !
    character(80)                         :: msgstr    ! string to pass message
    !
    logical                               :: cellfound ! indicate whether cell containing given point is found or not
    logical                               :: EQREAL    ! indicate whether two reals are equal or not
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanInterpolatePoint')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! assign exception value to output scalar (possibly overwritten by interpolated value)
    !
    foutp = excval
    !
    ! find closest vertex for given point
    !
    call SwanFindPoint ( x, y, ivert )
    !
    ! if point not found, give warning and return
    !
    if ( ivert < 0 ) then
       write (msgstr, '(a,f12.4,a,f12.4,a)') ' Point (',x+XOFFS,',',y+YOFFS,') not given in computational grid'
       call msgerr( 1, trim(msgstr) )
       return
    endif
    !
    ! if exception value found in closest vertex, return
    !
    if ( EQREAL(finp(ivert),excval) ) return
    !
    ! determine direction of given point to closest vertex
    !
    dxp = xcugrd(ivert) - x
    dyp = ycugrd(ivert) - y
    !
    ! if given point equals closest vertex, determine output quantity and return
    !
    if ( EQREAL(dxp,0.) .and. EQREAL(dyp,0.) ) then
       foutp = finp(ivert)
       return
    endif
    !
    th = atan2(dyp,dxp)
    !
    cellfound = .false.
    !
    ! loop over cells around closest vertex
    !
    celloop: do jc = 1, vert(ivert)%noc
       !
       ! get cell and its vertices
       !
       icell = vert(ivert)%cell(jc)%atti(CELLID)
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! get directions of faces to closest vertex
       !
       do k = 1, 3
          if ( v(k) == ivert ) then
             th1 = cell(icell)%geom(k)%th1
             th2 = cell(icell)%geom(k)%th2
             exit
          endif
       enddo
       !
       thdiff = th - th2
       do
          if ( abs(thdiff) <= PI ) exit
          th = th - sign (2., thdiff) * PI
          thdiff = th - th2
       enddo
       !
       ! is given point inside considered cell?
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) then   ! boundary vertex
          eps = PI/360.
       else
          eps = 0.
       endif
       !
       if ( th > th1-eps .and. th <= th2+eps ) then
          cellfound = .true.
          exit celloop
       endif
       !
    enddo celloop
    !
    ! if cell containing given point not found, give warning and return
    !
    if ( .not.cellfound ) then
       write (msgstr, '(a,f12.4,a,f12.4,a)') ' No triangle containing point (',x+XOFFS,',',y+YOFFS,') is found'
       call msgerr( 1, trim(msgstr) )
       return
    endif
    !
    ! determine output scalar in vertices
    !
    phi1 = finp(v(1))
    phi2 = finp(v(2))
    phi3 = finp(v(3))
    !
    ! 2D linear interpolation on considered triangle is carried out only if there are no exception values
    !
    if ( .not.EQREAL(phi1,excval) .and. .not.EQREAL(phi2,excval) .and. .not.EQREAL(phi3,excval) ) then
       !
       ! determine centroid and area of found cell
       !
       xc    = cell(icell)%attr(CELLCX  )
       yc    = cell(icell)%attr(CELLCY  )
       carea = cell(icell)%attr(CELLAREA)
       !
       ! determine output scalar in centroid
       !
       phic = ( phi1 + phi2 + phi3 ) / 3.
       !
       ! determine translation vectors of found cell
       !
       vec12(1) = xcugrd(v(2)) - xcugrd(v(1))
       vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
       vec23(1) = xcugrd(v(3)) - xcugrd(v(2))
       vec23(2) = ycugrd(v(3)) - ycugrd(v(2))
       vec31(1) = xcugrd(v(1)) - xcugrd(v(3))
       vec31(2) = ycugrd(v(1)) - ycugrd(v(3))
       !
       ! determine gradient scalar vector inside found cell based on outward normals
       ! Note: the outward normal is obtained by rotating the translation vector
       !       over 90 degrees in clockwise direction
       !
       xgrs =  vec23(2)*phi1 + vec31(2)*phi2 + vec12(2)*phi3
       ygrs = -vec23(1)*phi1 - vec31(1)*phi2 - vec12(1)*phi3
       !
       xgrs = -0.5*xgrs/carea
       ygrs = -0.5*ygrs/carea
       !
       ! determine output scalar inside considered triangle by means of 2D interpolation
       ! using constant gradient scalar vector
       !
       foutp = phic + xgrs*(x - xc) + ygrs*(y - yc)
       !
    endif
    !
end subroutine SwanInterpolatePoint
