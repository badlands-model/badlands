subroutine SwanTranspX ( amat   , rhs  , ac2   , ac1   , cax   , cay   , &
                         rdx    , rdy  , obredf, idcmin, idcmax, isslow, &
                         isstop , trac0, trac1 )
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
!   40.85, August 2008: add xy-propagation for output purposes
!
!   Purpose
!
!   Computes transport in x-y space using the lowest order upwind scheme
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
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
    real, dimension(MDC,MSC,nverts), intent(in) :: ac1    ! action density at previous time level
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(MDC,MSC,5), intent(inout)   :: amat   ! coefficient matrix of system of equations in (sigma,theta) space
                                                          ! 1: correspond to point (l  ,m  )
                                                          ! 2: correspond to point (l-1,m  )
                                                          ! 3: correspond to point (l+1,m  )
                                                          ! 4: correspond to point (l  ,m-1)
                                                          ! 5: correspond to point (l  ,m+1)
    real, dimension(MDC,MSC,ICMAX), intent(in)  :: cax    ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(in)  :: cay    ! wave transport velocity in y-direction
    real, dimension(MDC,MSC,2), intent(in)      :: obredf ! action reduction coefficient based on transmission
    real, dimension(2), intent(in)              :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)              :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MDC,MSC), intent(inout)     :: rhs    ! right-hand side of system of equations in (sigma,theta) space
    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac0  ! explicit part of propagation in present vertex for output purposes
    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac1  ! implicit part of propagation in present vertex for output purposes
!
!   Local variables
!
    integer       :: id       ! loop counter over direction bins
    integer       :: iddum    ! counter in directional space for considered sweep
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: is       ! loop counter over frequency bins
    integer       :: iv1      ! first index in computational stencil
    integer       :: iv2      ! second index in computational stencil
    integer       :: iv3      ! third index in computational stencil
    !
    real          :: acold    ! action density at previous time level
    real          :: asum     ! contributions to the matrix
    real          :: fac1     ! auxiliary factor
    real          :: fac2     ! another auxiliary factor
    real          :: rsum     ! contributions to the right-hand side
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanTranspX')
    !
    iv1 = vs(1)
    iv2 = vs(2)
    iv3 = vs(3)
    !
    if ( KSPHER == 0 ) then
       !
       fac1 = 1.
       fac2 = 1.
       !
    else
       !
       fac1 = COSLAT(2)/COSLAT(1)
       fac2 = COSLAT(3)/COSLAT(1)
       !
    endif
    !
    do is = isslow, isstop
       !
       do iddum = idcmin(is), idcmax(is)
          id = mod ( iddum - 1 + MDC , MDC ) + 1
          !
          ! compute the contributions based on the lowest order upwind scheme
          !
          asum = (rdx(1)+rdx(2))*cax(id,is,1) + (rdy(1)+rdy(2))*cay(id,is,1)
          !
          rsum = (rdx(1)*cax(id,is,2) + rdy(1)*cay(id,is,2)*fac1)*obredf(id,is,1)*ac2(id,is,iv2) + &
                 (rdx(2)*cax(id,is,3) + rdy(2)*cay(id,is,3)*fac2)*obredf(id,is,2)*ac2(id,is,iv3)
          !
          ! build the system
          !
          amat(id,is,1) = amat(id,is,1) + asum
          rhs (id,is  ) = rhs (id,is  ) + rsum
          !
          trac0(id,is,1) = trac0(id,is,1) - rsum
          trac1(id,is,1) = trac1(id,is,1) + asum
          !
          if ( NSTATC == 1 ) then
             !
             if ( ITERMX == 1 ) then
                acold = ac2(id,is,iv1)
             else
                acold = ac1(id,is,iv1)
             endif
             !
             amat(id,is,1) = amat(id,is,1) + RDTIM
             rhs (id,is  ) = rhs (id,is  ) + acold*RDTIM
             !
             trac0(id,is,1) = trac0(id,is,1) - acold*RDTIM
             trac1(id,is,1) = trac1(id,is,1) + RDTIM
             !
          endif
          !
       enddo
       !
    enddo
    !
end subroutine SwanTranspX
