subroutine SwanInterpolateOutput ( foutp, x, y, finp, mip, kvert, excval )
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
!   40.90: Nico Booij
!   41.07: Marcel Zijlema
!
!   Updates
!
!   40.80, August 2007: New subroutine
!   40.90,   June 2008: improved interpolation near obstacles
!   41.07,   July 2009: optimization
!
!   Purpose
!
!   Interpolates given output quantity to given points
!
!   Method
!
!   Look for closest vertex and determine triangle in which given point is located
!   Determine weighting coefficients for the corresponding vertices
!   Set weighting coeff to zero if there is an obstacle between given point and vertex
!   Interpolate output quantity using the resulting weighting coefficients
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use m_obsta
    use outp_data
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                 :: mip    ! number of given points
    integer, dimension(mip), intent(in) :: kvert  ! vertex indices of output points
    !
    real, intent(in)                    :: excval ! exception value for output quantity
    real, dimension(nverts), intent(in) :: finp   ! output quantity defined on the computational grid
    real, dimension(mip), intent(out)   :: foutp  ! interpolated output quantity at given points
    real, dimension(mip), intent(in)    :: x      ! x-coordinate of given points
    real, dimension(mip), intent(in)    :: y      ! y-coordinate of given points
!
!   Local variables
!
    integer                               :: icell     ! cell index
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: ip        ! loop counter
    integer, dimension(3)                 :: ivc       ! vertex indices in cyclic order
    integer                               :: ivert     ! vertex index
    integer                               :: jc        ! loop counter
    integer                               :: k         ! loop counter
    integer                               :: l         ! loop counter
    integer                               :: numcor    ! number of corner points in an obstacle
    integer, dimension(3)                 :: v         ! vertices in present cell
    !
    real                                  :: dxp       ! distance between given point and present vertex in x-direction
    real, dimension (3)                   :: dxv       ! difference of vertices of opposite side in x-coordinate
    real                                  :: dyp       ! distance between given point and present vertex in y-direction
    real, dimension (3)                   :: dyv       ! difference of vertices of opposite side in y-coordinate
    real                                  :: eps       ! a small number
    real                                  :: sumww     ! sum of the interpolation weights
    real                                  :: th        ! direction of given point to present vertex
    real                                  :: th1       ! direction of one face pointing to present vertex
    real                                  :: th2       ! direction of another face pointing to present vertex
    real                                  :: thdiff    ! difference between th and th2
    real                                  :: xb        ! user x-coordinate of begin of obstacle side
    real                                  :: xe        ! user x-coordinate of end of obstacle side
    real, dimension (3)                   :: xv        ! x-coordinate of the vertex
    real                                  :: yb        ! user y-coordinate of begin of obstacle side
    real                                  :: ye        ! user y-coordinate of end of obstacle side
    real, dimension (3)                   :: yv        ! y-coordinate of the vertex
    real, dimension (3)                   :: ww        ! weight of each vertex in the interpolation
    !
    character(80)                         :: msgstr    ! string to pass message
    !
    logical                               :: cellfound ! indicate whether cell containing given point is found or not
    logical, dimension (3)                :: cross     ! if true there is an obstacle between given point and vertex
    logical                               :: EQREAL    ! indicate whether two reals are equal or not
    logical                               :: obstcell  ! if true there is an obstacle in cell
    logical                               :: TCROSS    ! determines whether two line segments cross
    logical                               :: xonobst   ! not used
    !
    type(OBSTDAT), pointer                :: COBST     ! pointer to obstacle data
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
    if (ltrace) call strace (ient,'SwanInterpolateOutput')
    !
    if ( LCOMPGRD .and. mip == nverts ) then
       foutp = finp
       return
    endif
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! assign exception value to output quantity (possibly overwritten by interpolated values)
    !
    foutp = excval
    !
    ! loop over all given points
    !
    pointloop: do ip = 1, mip
       !
       ! assign vertex index of given point
       !
       ivert = kvert(ip)
       !
       ! if point not found, go to next point
       !
       if ( ivert < 0 ) cycle pointloop
       !
       ! if closest vertex is not active, go to next point
       !
       if ( .not.vert(ivert)%active ) cycle pointloop
       !
       ! determine direction of given point to closest vertex
       !
       dxp = xcugrd(ivert) - x(ip)
       dyp = ycugrd(ivert) - y(ip)
       !
       ! if given point equals closest vertex, determine output quantity and go to next point
       !
       if ( EQREAL(dxp,0.) .and. EQREAL(dyp,0.) ) then
          foutp(ip) = finp(ivert)
          cycle pointloop
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
       ! if cell containing given point not found, give warning and go to next point
       !
       if ( .not.cellfound ) then
          write (msgstr, '(a,f12.4,a,f12.4,a)') ' No triangle containing point (',x(ip)+XOFFS,',',y(ip)+YOFFS,') is found'
          call msgerr( 1, trim(msgstr) )
          cycle pointloop
       endif
       !
       ! 2D linear interpolation on considered triangle is carried out only if all vertices are active
       !
       if ( vert(v(1))%active .and. vert(v(2))%active .and. vert(v(3))%active ) then
          !
          !  get coordinates of the vertices
          !
          do k = 1, 3
             xv(k) = xcugrd(v(k))
             yv(k) = ycugrd(v(k))
             cross(k) = .false.
          enddo
          !
          ! determine difference in x and y of opposite side
          !
          do k = 1, 3
             ivc(2) = mod(k  ,3)+1
             ivc(3) = mod(k+1,3)+1
             dxv(k) = xv(ivc(3)) - xv(ivc(2))
             dyv(k) = yv(ivc(3)) - yv(ivc(2))
          enddo
          !
          ! determine whether there is an obstacle between given point and vertices
          !
          if ( NUMOBS > 0 ) then
             !
             COBST => FOBSTAC
             !
             do jc = 1, NUMOBS
                !
                numcor = COBST%NCRPTS
                if ( ITEST >= 120 ) write (PRINTF,10) jc, numcor
                !
                xb = COBST%XCRP(1)
                yb = COBST%YCRP(1)
                if ( ITEST >= 120 ) write (PRINTF,20) 1, xb+XOFFS, yb+YOFFS
                !
                do l = 2, numcor
                   !
                   xe = COBST%XCRP(l)
                   ye = COBST%YCRP(l)
                   if ( ITEST >= 120 ) write (PRINTF,20) l, xe+XOFFS, ye+YOFFS
                   !
                   ! loop over vertices
                   !
                   do k = 1, 3
                      if ( TCROSS(x(ip), xv(k), xb, xe, y(ip), yv(k), yb, ye, xonobst) ) cross(k) = .true.
                   enddo
                   !
                   xb = xe
                   yb = ye
                   !
                enddo
                !
                if (.not.associated(COBST%NEXTOBST)) exit
                COBST => COBST%NEXTOBST
                !
             enddo
             !
          endif
          !
          ! determine weighting coefficients
          !
          obstcell = .false.
          do k = 1, 3
             if (cross(k)) then
                ww(k) = 0.
                obstcell = .true.
             else
                ivc(1) = k
                ivc(2) = mod(k  ,3)+1
                ivc(3) = mod(k+1,3)+1
                ww(k) = ((x(ip) - xv(ivc(3))) * dyv(ivc(1)) - (y(ip) - yv(ivc(3))) * dxv(ivc(1))) / ( dxv(ivc(2)) * dyv(ivc(1)) - dyv(ivc(2)) * dxv(ivc(1)) )
             endif
          enddo
          if (obstcell) sumww = sum(ww)
          !
          ! use weighting coefficients to determine interpolated output quantity
          !
          do k = 1, 3
             if ( ww(k) > 1.e-10 ) then
                if (obstcell) ww(k) = ww(k) / sumww
             endif
          enddo
          foutp(ip) = sum (ww(:) * finp(v(:)))
          !
       endif
       !
    enddo pointloop
    !
 10 format (' Obstacle number : ', i4,'  has ', i4, ' corners')
 20 format (' Corner number:', i4,'    Xp: ', e10.4, ' Yp: ', e11.4)
    !
end subroutine SwanInterpolateOutput
