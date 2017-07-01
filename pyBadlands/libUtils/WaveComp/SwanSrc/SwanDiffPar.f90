subroutine SwanDiffPar ( ac2, dep2, spcsig )
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
!   41.02: Marcel Zijlema
!
!   Updates
!
!   41.02, February 2009: New subroutine
!
!   Purpose
!
!   Computes diffraction parameter and its derivatives in vertices
!
!   Method
!
!   Diffraction is approximated using the eikonal equation which
!   relates the wavenumber K to the separation factor k. Several
!   expressions of the eikonal equation have been presented in
!   the literature:
!                                 DH
!   Battjes (1968):   K^2 = k^2 + --
!                                 H
!
!   where H is the wave height and D is the Laplacian operator.
!
!                                 D.(pDH)
!   Berkhoff (1972):  K^2 = k^2 + -------
!                                   pH
!
!   where p = cc_g and D is the gradient operator in this case.
!
!   In both cases, the eikonal equation may be written as follows:
!
!   K = k (1+delta)^0.5
!
!   with
!
!           D.(pDH)
!   delta = -------
!           k^2 pH
!
!   From implementation point of view, the Battjes' eikonal
!   equation can be obtained if
!                                c_g = k
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use m_diffr
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2      ! action density at current time level
    real, dimension(nverts), intent(in)         :: dep2     ! water depth at current time level
    real, dimension(MSC), intent(in)            :: spcsig   ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell          ! index of present cell
    integer, save                         :: ient = 0       ! number of entries in this subroutine
    integer                               :: ivert          ! loop counter over vertices
    integer                               :: jc             ! loop counter
    integer                               :: jcell          ! index of next cell
    integer, parameter                    :: jeiko=1        ! choice parameter:
                                                            ! 0 = eikonal equation according to Battjes (1968)
                                                            ! 1 = eikonal equation according to Berkhoff (1972)
    !
    integer, dimension(3)                 :: v              ! vertices in present cell
    !
    real                                  :: area0          ! area of present cell
    real                                  :: area1          ! area of next cell
    double precision                      :: carea          ! twices the area of centroid dual around present vertex
    real                                  :: cg0            ! mean group velocity in centroid of present cell
    real                                  :: cg1            ! mean group velocity in centroid of next cell
    real                                  :: cslat          ! cosine of latitude
    real                                  :: ctot           ! zeroth moment of energy times group velocity
    real                                  :: delta          ! local diffraction parameter
    real                                  :: denom          ! a denominator
    real                                  :: deploc         ! local depth
    real                                  :: dgx0           ! x-component of diffusion gradient inside present cell
    real                                  :: dgx1           ! x-component of diffusion gradient inside next cell
    real                                  :: dgxdx          ! x-gradient of x-diffusion gradient component
    real                                  :: dgy0           ! y-component of diffusion gradient inside present cell
    real                                  :: dgy1           ! y-component of diffusion gradient inside next cell
    real                                  :: dgydy          ! y-gradient of y-diffusion gradient component
    real                                  :: dhsdx          ! x-gradient of wave height
    real                                  :: dhsdy          ! y-gradient of wave height
    real                                  :: difp0          ! diffraction parameter in centroid of present cell
    real                                  :: difp1          ! diffraction parameter in centroid of next cell
    real                                  :: etot           ! zeroth moment of the variance spectrum
    real                                  :: fmax           ! upper bound of frequency space for integration
    real                                  :: fmin           ! lower bound of frequency space for integration
    real                                  :: k0             ! mean wave number in centroid of present cell
    real                                  :: k1             ! mean wave number in centroid of next cell
    real                                  :: ktot           ! zeroth moment of energy times wave number
    double precision                      :: x0             ! x-coordinate of the centroid of present cell
    double precision                      :: x1             ! x-coordinate of the centroid of next cell
    double precision                      :: y0             ! y-coordinate of the centroid of present cell
    double precision                      :: y1             ! y-coordinate of the centroid of next cell
    !
    real, dimension(MSC)                  :: cgloc          ! group velocity
    real, dimension(MSC)                  :: kloc           ! wave number
    real, dimension(MSC)                  :: n              ! ratio of group and phase velocity
    real, dimension(MSC)                  :: nd             ! derivative of n with respect to depth
    real, dimension(MDC)                  :: ecs            ! help array containing (co)sine of spectral directions
    !
    real, dimension(:), allocatable       :: cg             ! mean group velocity
    real, dimension(:), allocatable       :: hs             ! wave height
    real, dimension(:), allocatable       :: k              ! mean wave number
    !
    real                                  :: SwanIntgratSpc ! integration of variance over a part of frequency space
    !
    type(celltype), dimension(:), pointer :: cell           ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert           ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanDiffPar')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! allocation and initialization of wave height and mean wave parameters
    !
    allocate(hs(nverts))
    allocate( k(nverts))
    allocate(cg(nverts))
    !
    hs = 0.
    k  = 10.
    cg = 0.
    !
    ! compute total energy, mean wave number and mean group velocity in vertices
    !
    do ivert = 1, nverts
       !
       deploc = dep2(ivert)
       !
       if ( deploc <= DEPMIN ) cycle
       !
       ! compute group velocity and wave number for all frequencies
       !
       call KSCIP1 (MSC,spcsig,deploc,kloc,cgloc,n,nd)
       !
       ! integration over f in [0,infty]
       !
       fmin = 0.
       fmax = 1000.
       ecs  = 1.
       !
       etot = SwanIntgratSpc(0.   , fmin, fmax, spcsig, ecs,            &
                             kloc , ecs , 0.  , 0.    , ac2(:,:,ivert), &
                             1    )
       !
       ktot = SwanIntgratSpc(1.   , fmin, fmax, spcsig, ecs,            &
                             kloc , ecs , 0.  , 0.    , ac2(:,:,ivert), &
                             3    )
       !
       ctot = SwanIntgratSpc(1.   , fmin, fmax, spcsig, ecs,            &
                             cgloc, ecs , 0.  , 0.    , ac2(:,:,ivert), &
                             4    )
       !
       if ( etot > 0. ) then
          hs(ivert) = 4.*sqrt(etot)
          k (ivert) = ktot/etot
          cg(ivert) = ctot/etot
       endif
       !
    enddo
    !
    if ( jeiko == 0 ) cg = k
    !
    ! compute diffraction parameter in vertices
    !
    DIFPARAM = 1.
    !
    vertexloop : do ivert = 1, nverts
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop    ! boundary vertex
       !
       cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
       !
       ! compute contributions to the Laplacian in present vertex
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
          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop
          !
          ! determine centroid of present cell
          !
          x0    = cell(icell)%attr(CELLCX)
          y0    = cell(icell)%attr(CELLCY)
          area0 = cell(icell)%attr(CELLAREA)
          !
          ! determine mean wave number and group velocity in centroid in present cell
          !
          cg0 = ( cg(v(1)) + cg(v(2)) + cg(v(3)) )/ 3.
          k0  = ( k (v(1)) + k (v(2)) + k (v(3)) )/ 3.
          !
          ! determine derivatives of wave height inside present cell
          !
          dhsdx = 0.5*( hs(v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
                        hs(v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
                        hs(v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area0
          !
          dhsdy = 0.5*( hs(v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
                        hs(v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
                        hs(v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area0
          !
          ! in case of spherical coordinates, transform back to Cartesian coordinates
          !
          if ( KSPHER > 0 ) then
             !
             dhsdx = dhsdx/(cslat * LENDEG)
             dhsdy = dhsdy/LENDEG
             !
          endif
          !
          ! determine diffusion gradients in centroid of present cell
          !
          dgx0 = cg0*dhsdx/k0
          dgy0 = cg0*dhsdy/k0
          !
          ! get next cell in counterclockwise direction
          !
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          !
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          !
          ! determine centroid of next cell
          !
          x1    = cell(jcell)%attr(CELLCX)
          y1    = cell(jcell)%attr(CELLCY)
          area1 = cell(jcell)%attr(CELLAREA)
          !
          ! determine mean wave number and group velocity in centroid in next cell
          !
          cg1 = ( cg(v(1)) + cg(v(2)) + cg(v(3)) )/ 3.
          k1  = ( k (v(1)) + k (v(2)) + k (v(3)) )/ 3.
          !
          ! determine derivatives of wave height inside next cell
          !
          dhsdx = 0.5*( hs(v(1))*(ycugrd(v(2))-ycugrd(v(3))) + &
                        hs(v(2))*(ycugrd(v(3))-ycugrd(v(1))) + &
                        hs(v(3))*(ycugrd(v(1))-ycugrd(v(2))) )/area1
          !
          dhsdy = 0.5*( hs(v(1))*(xcugrd(v(3))-xcugrd(v(2))) + &
                        hs(v(2))*(xcugrd(v(1))-xcugrd(v(3))) + &
                        hs(v(3))*(xcugrd(v(2))-xcugrd(v(1))) )/area1
          !
          ! in case of spherical coordinates, transform back to Cartesian coordinates
          !
          if ( KSPHER > 0 ) then
             !
             dhsdx = dhsdx/(cslat * LENDEG)
             dhsdy = dhsdy/LENDEG
             !
          endif
          !
          ! determine diffusion gradients in centroid of next cell
          !
          dgx1 = cg1*dhsdx/k1
          dgy1 = cg1*dhsdy/k1
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
          denom = k(ivert)*cg(ivert)*hs(ivert)
          !
          if ( denom > 0. ) then
             delta = (dgxdx + dgydy)/denom
          else
             delta = 0.
          endif
          !
          if ( delta > -1. ) DIFPARAM(ivert) = sqrt(1.+delta)
          !
       endif
       !
    enddo vertexloop
    !
    ! deallocation of wave parameters
    !
    deallocate(cg)
    deallocate(hs)
    deallocate(k )
    !
    ! compute derivatives of diffraction parameter in vertices
    !
    DIFPARDX = 0.
    DIFPARDY = 0.
    !
    vertexloop2 : do ivert = 1, nverts
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop2   ! boundary vertex
       !
       cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
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
          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop2
          !
          ! determine centroid of present cell
          !
          x0 = cell(icell)%attr(CELLCX)
          y0 = cell(icell)%attr(CELLCY)
          !
          ! determine diffraction parameter in centroid in present cell
          !
          difp0 = ( DIFPARAM(v(1)) + DIFPARAM(v(2)) + DIFPARAM(v(3)) )/ 3.
          !
          ! get next cell in counterclockwise direction
          !
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          !
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          !
          ! determine centroid of next cell
          !
          x1 = cell(jcell)%attr(CELLCX)
          y1 = cell(jcell)%attr(CELLCY)
          !
          ! determine diffraction parameter in centroid of next cell
          !
          difp1 = ( DIFPARAM(v(1)) + DIFPARAM(v(2)) + DIFPARAM(v(3)) )/ 3.
          !
          ! compute contribution to area of centroid dual
          !
          carea = carea + x0*y1 - x1*y0
          !
          ! compute x-gradient of diffraction parameter
          !
          dgxdx = dgxdx + ( difp0 + difp1 ) * real( y1 - y0 )
          !
          ! compute y-gradient of diffraction parameter
          !
          dgydy = dgydy + ( difp0 + difp1 ) * real( x0 - x1 )
          !
       enddo
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
          DIFPARDX(ivert) = dgxdx
          DIFPARDY(ivert) = dgydy
          !
       endif
       !
    enddo vertexloop2
    !
end subroutine SwanDiffPar
