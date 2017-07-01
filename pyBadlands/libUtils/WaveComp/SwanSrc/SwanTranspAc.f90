subroutine SwanTranspAc ( amat  , rhs   , leakcf, ac2   , ac1   , &
                          cgo   , cax   , cay   , cad   , cas   , &
                          anybin, rdx   , rdy   , spcsig, spcdir, &
                          obredf, idcmin, idcmax, iscmin, iscmax, &
                          iddlow, iddtop, isslow, isstop, anyblk, &
                          trac0 , trac1 )
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
!   40.80,   August 2007: New subroutine
!   40.85,   August 2008: add tranport terms for output purposes
!   41.00, February 2009: add GSE correction
!   41.07,     July 2009: add explicit scheme for sigma space
!
!   Purpose
!
!   Computes the transport part of the action balance equation
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                         :: iddlow ! minimum direction bin that is propagated within a sweep
    integer, intent(in)                         :: iddtop ! maximum direction bin that is propagated within a sweep
    integer, intent(in)                         :: isslow ! minimum frequency that is propagated within a sweep
    integer, intent(in)                         :: isstop ! maximum frequency that is propagated within a sweep
    !
    integer, dimension(MSC), intent(in)         :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(in)         :: idcmin ! minimum frequency-dependent counter in directional space
    integer, dimension(MDC), intent(in)         :: iscmax ! maximum direction-dependent counter in frequency space
    integer, dimension(MDC), intent(in)         :: iscmin ! minimum direction-dependent counter in frequency space
    !
    real, dimension(MDC,MSC,nverts), intent(in) :: ac1    ! action density at previous time level
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(MDC,MSC,5), intent(out)     :: amat   ! coefficient matrix of system of equations in (sigma,theta) space
                                                          ! 1: correspond to point (l  ,m  )
                                                          ! 2: correspond to point (l-1,m  )
                                                          ! 3: correspond to point (l+1,m  )
                                                          ! 4: correspond to point (l  ,m-1)
                                                          ! 5: correspond to point (l  ,m+1)
    real, dimension(MDC,MSC), intent(in)        :: cad    ! wave transport velocity in theta-direction
    real, dimension(MDC,MSC), intent(in)        :: cas    ! wave transport velocity in sigma-direction
    real, dimension(MDC,MSC,ICMAX), intent(in)  :: cax    ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(in)  :: cay    ! wave transport velocity in y-direction
    real, dimension(MSC,ICMAX), intent(in)      :: cgo    ! group velocity
    real, dimension(MDC,MSC), intent(out)       :: leakcf ! leak coefficient
    real, dimension(MDC,MSC,2), intent(in)      :: obredf ! action reduction coefficient based on transmission
    real, dimension(2), intent(in)              :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)              :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MDC,MSC), intent(out)       :: rhs    ! right-hand side of system of equations in (sigma,theta) space
    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
                                                          ! (*,2): cosine of spectral directions
                                                          ! (*,3): sine of spectral directions
                                                          ! (*,4): cosine^2 of spectral directions
                                                          ! (*,5): cosine*sine of spectral directions
                                                          ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac0  ! explicit part of propagation in present vertex for output purposes
    real, dimension(MDC,MSC,MTRNP), intent(out) :: trac1  ! implicit part of propagation in present vertex for output purposes
    !
    logical, dimension(MDC,MSC), intent(in)     :: anybin ! true if bin is active in considered sweep
    logical, dimension(MDC,MSC), intent(out)    :: anyblk ! true if bin is blocked by a counter current based on a CFL criterion
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanTranspAc')
    !
    ! set matrix and right-hand side to zero
    !
    amat = 0.
    rhs  = 0.
    !
    ! set leak and transport coefficients to zero
    !
    leakcf = 0.
    trac0  = 0.
    trac1  = 0.
    !
    ! compute transport in x-y space
    !
!TIMG    call SWTSTA(140)
    call SwanTranspX ( amat   , rhs  , ac2   , ac1   , cax   , cay   , &
                       rdx    , rdy  , obredf, idcmin, idcmax, isslow, &
                       isstop , trac0, trac1 )
    !
    ! add GSE correction, if appropriate
    !
    if ( WAVAGE > 0. ) call SwanGSECorr ( rhs, ac2, cgo, spcdir, idcmin, idcmax, isslow, isstop, trac0 )
!TIMG    call SWTSTO(140)
    !
    ! compute transport in theta space
    !
!TIMG    call SWTSTA(142)
    if ( IREFR /= 0 ) then
       !
       call STRSD ( DDIR       , idcmin     , idcmax     , cad    , &
                    amat(1,1,4), amat(1,1,1), amat(1,1,5), rhs    , &
                    ac2        , isstop     , anybin     , leakcf , &
                    trac0      , trac1      )
       !
    endif
!TIMG    call SWTSTO(142)
    !
    ! compute transport in sigma space
    !
!TIMG    call SWTSTA(141)
    if ( (DYNDEP .OR. ICUR /= 0) .and. ITFRE /= 0 ) then
       !
       if ( int(PNUMS(8)) == 1 ) then
          !
          ! implicit scheme
          !
          call STRSSI ( spcsig     , cas   , amat(1,1,2), amat(1,1,1), &
                        amat(1,1,3), anybin, rhs        , ac2        , &
                        iscmin     , iscmax, iddlow     , iddtop     , &
                        trac0      , trac1 )
          !
       elseif ( int(PNUMS(8)) == 2 ) then
          !
          ! explicit scheme
          !
          call STRSSB ( iddlow, iddtop, idcmin, idcmax, isstop, &
                        cax   , cay   , cas   , ac2   , spcsig, &
                        rhs   , anyblk, rdx   , rdy   , trac0 )
          !
       endif
       !
    endif
!TIMG    call SWTSTO(141)
    !
end subroutine SwanTranspAc
