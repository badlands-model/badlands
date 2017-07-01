subroutine SwanComputeForce ( fx, fy, ac2, dep2, spcsig, spcdir )
!ADCsubroutine SwanComputeForce
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
!   41.20: Casey Dietrich
!
!   Updates
!
!   40.80, February 2008: New subroutine
!   41.20, February 2008: extension to tightly coupled ADCIRC+SWAN model
!
!   Purpose
!
!   Computes wave-induced force in vertices
!
!   Method
!
!   First, compute radiation stresses in vertices
!   Next, compute gradients of radiation stresses in vertices
!   Finally, compute wave-induced force in vertices
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects
!ADC    use Couple2Adcirc, only: compda
!ADC    use GLOBAL,        only: rsnx2, rsny2
!ADC    use m_genarr,      only: ac2, spcdir, spcsig
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(nverts), intent(in)         :: dep2   ! water depth at current time level
    real, dimension(nverts), intent(out)        :: fx     ! wave-induced force in x-direction
    real, dimension(nverts), intent(out)        :: fy     ! wave-induced force in y-direction
    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
                                                          ! (*,2): cosine of spectral directions
                                                          ! (*,3): sine of spectral directions
                                                          ! (*,4): cosine^2 of spectral directions
                                                          ! (*,5): cosine*sine of spectral directions
                                                          ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell    ! index of present cell
    integer                               :: id       ! loop counter over direction bins
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: jc       ! loop counter
    integer                               :: jcell    ! index of next cell
    !
    integer, dimension(3)                 :: v        ! vertices in present cell
    !
    double precision                      :: area     ! twices the area of centroid dual around present vertex
    real                                  :: cslat    ! cosine of latitude
    real                                  :: deploc   ! local depth
    real                                  :: dsxxdx   ! x-gradient of sxx
    real                                  :: dsxydx   ! x-gradient of sxy
    real                                  :: dsxydy   ! y-gradient of sxy
    real                                  :: dsyydy   ! y-gradient of syy
    real                                  :: sxx0     ! sxx in centroid of present cell
    real                                  :: sxx1     ! sxx in centroid of next cell
    real                                  :: sxy0     ! sxy in centroid of present cell
    real                                  :: sxy1     ! sxy in centroid of next cell
    real                                  :: syy0     ! syy in centroid of present cell
    real                                  :: syy1     ! syy in centroid of next cell
    real                                  :: sxxsumd  ! cumulated sxx over direction space
    real                                  :: sxysumd  ! cumulated sxy over direction space
    real                                  :: syysumd  ! cumulated syy over direction space
    real                                  :: sxxsums  ! cumulated sxx over frequency space
    real                                  :: sxysums  ! cumulated sxy over frequency space
    real                                  :: syysums  ! cumulated syy over frequency space
    double precision                      :: x0       ! x-coordinate of the centroid of present cell
    double precision                      :: x1       ! x-coordinate of the centroid of next cell
    double precision                      :: y0       ! y-coordinate of the centroid of present cell
    double precision                      :: y1       ! y-coordinate of the centroid of next cell
    !
    real, dimension(1)                    :: cg       ! group velocity
    real, dimension(1)                    :: k        ! wave number
    real, dimension(1)                    :: n        ! ratio of group and phase velocity
    real, dimension(1)                    :: nd       ! derivative of n with respect to depth
    real, dimension(1)                    :: sig      ! relative frequency
!ADC    !
!ADC    real, dimension(nverts)               :: dep2     ! help array to store water depth
!ADC    real, dimension(nverts)               :: fx       ! help array to store wave-induced force in x-direction
!ADC    real, dimension(nverts)               :: fy       ! help array to store wave-induced force in y-direction
    !
    real, dimension(:), allocatable       :: sxx      ! x-component of radiation stress in x-direction
    real, dimension(:), allocatable       :: sxy      ! cross component of radiation stress in x/y-direction
    real, dimension(:), allocatable       :: syy      ! y-component of radiation stress in y-direction
    !
    character(80)                         :: msgstr   ! string to pass message
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
    if (ltrace) call strace (ient,'SwanComputeForce')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! allocation and initialization of radiation stresses
    !
    allocate(sxx(nverts))
    allocate(sxy(nverts))
    allocate(syy(nverts))
    !
    sxx = 0.
    sxy = 0.
    syy = 0.
!ADC    !
!ADC    ! move the depths to their own array
!ADC    do ivert = 1, nverts
!ADC       dep2(ivert) = compda(ivert,jdp2)
!ADC    enddo
    !
    ! compute radiation stresses in vertices
    !
    do ivert = 1, nverts
       !
       deploc = dep2(ivert)
       !
       if ( deploc <= DEPMIN ) cycle
       !
       sxxsums = 0.
       sxysums = 0.
       syysums = 0.
       !
       ! calculate sum of each component contributed to the radiation stresses over all frequencies
       !
       do is = 1, MSC
          !
          ! compute group velocity over phase velocity (=n)
          !
          sig(1) = spcsig(is)
          call KSCIP1 (1,sig,deploc,k,cg,n,nd)
          !
          sxxsumd = 0.
          sxysumd = 0.
          syysumd = 0.
          !
          ! calculate sum of each component contributed to the radiation stresses over all wave directions
          !
          do id = 1, MDC
             !
             sxxsumd = sxxsumd + ( n(1)*(spcdir(id,4) + 1.) -0.5 ) * ac2(id,is,ivert)
             sxysumd = sxysumd +   n(1)* spcdir(id,5)              * ac2(id,is,ivert)
             syysumd = syysumd + ( n(1)*(spcdir(id,6) + 1.) -0.5 ) * ac2(id,is,ivert)
             !
          enddo
          !
          ! integrate over frequency space
          !
          sxxsums = sxxsums + sxxsumd * spcsig(is)**2
          sxysums = sxysums + sxysumd * spcsig(is)**2
          syysums = syysums + syysumd * spcsig(is)**2
          !
       enddo
       !
       ! compute the radiation stresses (divided by rho times gravitational acceleration)
       !
       sxx(ivert) = sxxsums * DDIR * FRINTF
       sxy(ivert) = sxysums * DDIR * FRINTF
       syy(ivert) = syysums * DDIR * FRINTF
       !
    enddo
    !
    ! compute wave-induced force in vertices
    !
    fx = 0.
    fy = 0.
    !
    vertexloop : do ivert = 1, nverts
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop    ! boundary vertex
       !
       ! first, compute gradients of the radiation stresses in vertices
       !
       area   = 0d0
       dsxxdx = 0.
       dsxydx = 0.
       dsxydy = 0.
       dsyydy = 0.
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
          x0 = cell(icell)%attr(CELLCX)
          y0 = cell(icell)%attr(CELLCY)
          !
          ! determine radiation stresses in centroid in present cell
          !
          sxx0 = ( sxx(v(1)) + sxx(v(2)) + sxx(v(3)) )/ 3.
          sxy0 = ( sxy(v(1)) + sxy(v(2)) + sxy(v(3)) )/ 3.
          syy0 = ( syy(v(1)) + syy(v(2)) + syy(v(3)) )/ 3.
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
          ! determine radiation stresses in centroid of next cell
          !
          sxx1 = ( sxx(v(1)) + sxx(v(2)) + sxx(v(3)) )/ 3.
          sxy1 = ( sxy(v(1)) + sxy(v(2)) + sxy(v(3)) )/ 3.
          syy1 = ( syy(v(1)) + syy(v(2)) + syy(v(3)) )/ 3.
          !
          ! compute contribution to area of centroid dual
          !
          area = area + x0*y1 - x1*y0
          !
          ! compute contribution to x-gradient of radiation stresses sxx and sxy
          !
          dsxxdx = dsxxdx + ( sxx0 + sxx1 ) * real( y1 - y0 )
          dsxydx = dsxydx + ( sxy0 + sxy1 ) * real( y1 - y0 )
          !
          ! compute contribution to y-gradient of radiation stresses sxy and syy
          !
          dsxydy = dsxydy + ( sxy0 + sxy1 ) * real( x1 - x0 )
          dsyydy = dsyydy + ( syy0 + syy1 ) * real( x1 - x0 )
          !
       enddo
       !
       ! if area is non-positive, give error and go to next vertex
       !
       if ( area <= 0d0 ) then
          write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
          call msgerr( 2, trim(msgstr) )
          cycle vertexloop
       endif
       !
       dsxxdx =  dsxxdx/real(area)
       dsxydx =  dsxydx/real(area)
       dsxydy = -dsxydy/real(area)
       dsyydy = -dsyydy/real(area)
       !
       ! in case of spherical coordinates, transform back to Cartesian coordinates
       !
       if ( KSPHER > 0 ) then
          !
          cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
          !
          dsxxdx = dsxxdx/(cslat * LENDEG)
          dsxydx = dsxydx/(cslat * LENDEG)
          dsxydy = dsxydy/LENDEG
          dsyydy = dsyydy/LENDEG
          !
       endif
       !
       ! finally, compute wave-induced force
       !
       fx(ivert) = -RHO * GRAV * ( dsxxdx  + dsxydy )
       fy(ivert) = -RHO * GRAV * ( dsxydx  + dsyydy )
       !
    enddo vertexloop
!ADC    !
!ADC    ! pass stresses to ADCIRC
!ADC    do ivert = 1, nverts
!ADC       rsnx2(ivert) = fx(ivert)
!ADC       rsny2(ivert) = fy(ivert)
!ADC    enddo
    !
    ! deallocation of radiation stresses
    !
    deallocate(sxx)
    deallocate(sxy)
    deallocate(syy)
    !
end subroutine SwanComputeForce
