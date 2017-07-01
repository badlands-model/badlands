subroutine SwanGSECorr ( rhs, ac2, cgo, spcdir, idcmin, idcmax, isslow, isstop, trac0 )
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
!   41.00: Marcel Zijlema
!
!   Updates
!
!   41.00, February 2009: New subroutine
!
!   Purpose
!
!   Computes waveage-dependent diffusion terms in x-y space to counteract the garden-sprinkler effect
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                         :: isslow ! minimum frequency that is propagated within a sweep
    integer, intent(in)                         :: isstop ! maximum frequency that is propagated within a sweep
    !
    integer, dimension(MSC), intent(in)         :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(in)         :: idcmin ! minimum frequency-dependent counter in directional space
    !
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(MSC,ICMAX), intent(in)      :: cgo    ! group velocity
    real, dimension(MDC,MSC), intent(inout)     :: rhs    ! right-hand side of system of equations in (sigma,theta) space
    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
                                                          ! (*,2): cosine of spectral directions
                                                          ! (*,3): sine of spectral directions
                                                          ! (*,4): cosine^2 of spectral directions
                                                          ! (*,5): cosine*sine of spectral directions
                                                          ! (*,6): sine^2 of spectral directions
    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac0  ! explicit part of propagation in present vertex for output purposes
!
!   Local variables
!
    integer                               :: icell    ! index of present cell
    integer                               :: id       ! loop counter over direction bins
    integer                               :: iddum    ! counter in directional space for considered sweep
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: ivert    ! index of present vertex
    integer                               :: jc       ! loop counter
    integer                               :: jcell    ! index of next cell
    !
    integer, dimension(3)                 :: v        ! vertices in present cell
    !
    real                                  :: area0    ! area of present cell
    real                                  :: area1    ! area of next cell
    double precision                      :: carea    ! twices the area of centroid dual around present vertex
    real                                  :: cslat    ! cosine of latitude
    real                                  :: dac2dx   ! x-derivative of action density
    real                                  :: dac2dy   ! y-derivative of action density
    real                                  :: dcg      ! group velocity difference across frequency bin
    real                                  :: dgx0     ! x-component of diffusion gradient inside present cell
    real                                  :: dgx1     ! x-component of diffusion gradient inside next cell
    real                                  :: dgxdx    ! x-gradient of x-diffusion gradient component
    real                                  :: dgy0     ! y-component of diffusion gradient inside present cell
    real                                  :: dgy1     ! y-component of diffusion gradient inside next cell
    real                                  :: dgydy    ! y-gradient of y-diffusion gradient component
    real                                  :: dnn      ! waveage-dependent diffusion coefficient normal to propagation direction
    real                                  :: dss      ! waveage-dependent diffusion coefficient in propagation direction
    real                                  :: dxx      ! xx-component of diffusion coefficient in Cartesian coordinates
    real                                  :: dxy      ! xy-component of diffusion coefficient in Cartesian coordinates
    real                                  :: dyy      ! yy-component of diffusion coefficient in Cartesian coordinates
    double precision                      :: x0       ! x-coordinate of the centroid of present cell
    double precision                      :: x1       ! x-coordinate of the centroid of next cell
    double precision                      :: y0       ! y-coordinate of the centroid of present cell
    double precision                      :: y1       ! y-coordinate of the centroid of next cell
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGSECorr')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ivert = vs(1)
    !
    if ( vert(ivert)%atti(VMARKER) == 1 ) return    ! no GSE correction in boundary vertex
    !
    cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
    !
    do is = isslow, isstop
       !
       ! calculate waveage-dependent diffusion coefficients in polar coordinates
       !
       if ( is == 1 ) then
          dcg = abs(cgo(is+1,1)-cgo(is,1))
       elseif ( is == isstop ) then
          dcg = abs(cgo(is,1)-cgo(is-1,1))
       else
          dcg = 0.5 * abs(cgo(is+1,1)-cgo(is-1,1))
       endif
       !
       dss = dcg**2*WAVAGE/12.
       dnn = (cgo(is,1)*DDIR)**2 * WAVAGE/12.
       !
       do iddum = idcmin(is), idcmax(is)
          id = mod ( iddum - 1 + MDC , MDC ) + 1
          !
          ! calculate diffusion coefficients in Cartesian coordinates
          !
          dxx = dss*spcdir(id,4) + dnn*spcdir(id,6)
          dyy = dss*spcdir(id,6) + dnn*spcdir(id,4)
          dxy = (dss-dnn)*spcdir(id,5)
          !
          ! compute contribution to the diffusion terms in present vertex
          !
          carea = 0d0
          dgxdx = 0.
          dgydy = 0.
          !
          ! loop over cells around considered vertex
          !
          do jc = 1, vert(ivert)%noc
             !
             ! get present cell and its vertices
             !
             icell = vert(ivert)%cell(jc)%atti(CELLID)
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             ! determine centroid and area of present cell
             !
             x0    = cell(icell)%attr(CELLCX)
             y0    = cell(icell)%attr(CELLCY)
             area0 = cell(icell)%attr(CELLAREA)
             !
             ! determine derivatives of action density inside present cell
             !
             dac2dx = 0.5*( ac2(id,is,v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
                            ac2(id,is,v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
                            ac2(id,is,v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area0
             !
             dac2dy = 0.5*( ac2(id,is,v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
                            ac2(id,is,v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
                            ac2(id,is,v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area0
             !
             ! in case of spherical coordinates, transform back to Cartesian coordinates
             !
             if ( KSPHER > 0 ) then
                !
                dac2dx = dac2dx/(cslat * LENDEG)
                dac2dy = dac2dy/LENDEG
                !
             endif
             !
             ! determine diffusion gradients in centroid of present cell
             !
             dgx0 = dxx*dac2dx + dxy*dac2dy
             dgy0 = dxy*dac2dx + dyy*dac2dy
             !
             ! get next cell in counterclockwise direction
             !
             jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
             !
             v(1) = cell(jcell)%atti(CELLV1)
             v(2) = cell(jcell)%atti(CELLV2)
             v(3) = cell(jcell)%atti(CELLV3)
             !
             ! determine centroid and area of next cell
             !
             x1    = cell(jcell)%attr(CELLCX)
             y1    = cell(jcell)%attr(CELLCY)
             area1 = cell(jcell)%attr(CELLAREA)
             !
             ! determine derivatives of action density inside next cell
             !
             dac2dx = 0.5*( ac2(id,is,v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
                            ac2(id,is,v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
                            ac2(id,is,v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area1
             !
             dac2dy = 0.5*( ac2(id,is,v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
                            ac2(id,is,v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
                            ac2(id,is,v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area1
             !
             ! in case of spherical coordinates, transform back to Cartesian coordinates
             !
             if ( KSPHER > 0 ) then
                !
                dac2dx = dac2dx/(cslat * LENDEG)
                dac2dy = dac2dy/LENDEG
                !
             endif
             !
             ! determine diffusion gradients in centroid of next cell
             !
             dgx1 = dxx*dac2dx + dxy*dac2dy
             dgy1 = dxy*dac2dx + dyy*dac2dy
             !
             ! compute contribution to area of centroid dual
             !
             carea = carea + x0*y1 - x1*y0
             !
             ! compute contribution to x-gradient of x-diffusion gradient in centroid dual
             !
             dgxdx = dgxdx + ( dgx0 + dgx1 ) * real( y1 - y0 )
             !
             ! compute contribution to y-gradient of y-diffusion gradient in centroid dual
             !
             dgydy = dgydy + ( dgy0 + dgy1 ) * real( x0 - x1 )
             !
          enddo
          !
          ! add diffusion terms to the right-hand side of the action balance equation
          !
          if ( carea > 0d0 ) then
             !
             dgxdx = dgxdx/real(carea)
             dgydy = dgydy/real(carea)
             !
             ! in case of spherical coordinates, transform back to Cartesian coordinates
             !
             if ( KSPHER > 0 ) then
                !
                dgxdx = dgxdx/(cslat * LENDEG)
                dgydy = dgydy/LENDEG
                !
             endif
             !
             rhs  (id,is  ) = rhs  (id,is  ) + dgxdx + dgydy
             trac0(id,is,1) = trac0(id,is,1) - dgxdx - dgydy
             !
          endif
          !
       enddo
       !
    enddo
    !
end subroutine SwanGSECorr
