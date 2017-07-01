subroutine SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, delhs, deltm, xytst, spcsig, ac2, ivlow, ivup )
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
!   40.93: Marcel Zijlema
!   41.10: Marcel Zijlema
!
!   Updates
!
!   40.80,   October 2007: New subroutine
!   40.93, September 2008: extended with curvature of mean period
!   41.10,    August 2009: parallelization using OpenMP directives
!
!   Purpose
!
!   Determine accuracy of wave height and period by means of curvature for convergence check
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                         :: ivlow  ! lower index in range of vertices in calling thread
    integer, intent(in)                         :: ivup   ! upper index in range of vertices in calling thread
    integer, dimension(NPTST), intent(in)       :: xytst  ! test points for output purposes
    !
    real, intent(out)                           :: accur  ! percentage of active vertices in which required accuracy has been reached
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(nverts), intent(out)        :: delhs  ! difference in wave height between last 2 iterations in all vertices
    real, dimension(nverts), intent(out)        :: deltm  ! difference in mean period between last 2 iterations in all vertices
    real, dimension(nverts), intent(inout)      :: hscurr ! wave height at current iteration level
    real, dimension(nverts), intent(inout)      :: hsdifc ! difference in wave height of current and one before previous iteration
    real, dimension(nverts), intent(inout)      :: hsprev ! wave height at previous iteration level
    real, dimension(nverts), intent(inout)      :: tmcurr ! mean period at current iteration level
    real, dimension(nverts), intent(inout)      :: tmdifc ! difference in mean period of current and one before previous iteration
    real, dimension(nverts), intent(inout)      :: tmprev ! mean period at previous iteration level
    real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: id       ! loop counter over direction bins
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: j        ! loop counter
    !
    real                                  :: curvah   ! required accuracy with respect to curvature in wave height
    real                                  :: curvat   ! required accuracy with respect to curvature in mean period
    real                                  :: fact     ! auxiliary factor
    real                                  :: hsabs    ! absolute difference in wave height between last 2 iterations
    real                                  :: hscurv   ! curvature of iteration curve of wave height
    real                                  :: hsdif0   ! value of hsdifc at previous iteration level
    real                                  :: hsprev0  ! wave height at one before previous iteration level
    real                                  :: hsrel    ! required accuracy with respect to relative error in wave height
    real                                  :: m0       ! moment of zeroth order
    real                                  :: m1       ! moment of first order
    real                                  :: npacc    ! number of vertices in which required accuracy has been reached
    real                                  :: npacct   ! npacc counter for calling thread
    real                                  :: nwetp    ! total number of active vertices
    real                                  :: nwetpt   ! nwetp counter for calling thread
    real                                  :: tmabs    ! absolute difference in mean period between last 2 iterations
    real                                  :: tmcurv   ! curvature of iteration curve of mean period
    real                                  :: tmdif0   ! value of tmdifc at previous iteration level
    real                                  :: tmprev0  ! mean period at one before previous iteration level
    real                                  :: tmrel    ! required accuracy with respect to relative error in mean period
    !
    logical                               :: lhead    ! logical indicating to write header
    logical                               :: tstfl    ! indicates whether vertex is a test point
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Common variables
!
    common/convstopc/npacc,nwetp
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanConvStopc')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    !$omp single
    npacc = 0.
    nwetp = 0.
    !$omp end single
    !
    npacct = 0.
    nwetpt = 0.
    !
    deltm = 0.
    delhs = 0.
    !
    lhead = .true.
    !
    ! calculate a set of accuracy parameters based on relative error and curvature for Hs and Tm
    !
    do ivert = ivlow, ivup
       !
       if ( vert(ivert)%active ) then
          !
          ! determine whether the present vertex is a test point
          !
          tstfl = .false.
          if ( NPTST > 0 ) then
             do j = 1, NPTST
                if ( ivert /= xytst(j) ) cycle
                tstfl = .true.
             enddo
          endif
          !
          ! count active points
          !
          nwetpt = nwetpt + 1.
          !
          ! store wave height and mean period of previous iteration levels
          !
          hsprev0       = max( 1.e-20, hsprev(ivert) )
          hsprev(ivert) = max( 1.e-20, hscurr(ivert) )
          tmprev0       = max( 1.e-20, tmprev(ivert) )
          tmprev(ivert) = max( 1.e-20, tmcurr(ivert) )
          !
          ! compute wave height and mean period for present vertex
          !
          m0 = 0.
          m1 = 0.
          do is = 1, MSC
             do id = 1, MDC
                fact = spcsig(is)**2 * ac2(id,is,ivert)
                m0 = m0 + fact
                m1 = m1 + fact * spcsig(is)
             enddo
          enddo
          m0 = m0 * FRINTF * DDIR
          m1 = m1 * FRINTF * DDIR
          !
          if ( m0 > 0. ) then
             hscurr(ivert) = max ( 1.e-20, 4.*sqrt(m0) )
          else
             hscurr(ivert) = 1.e-20
          endif
          if ( m1 > 0. ) then
             tmcurr(ivert) = max ( 1.e-20, PI2*(m0/m1) )
          else
             tmcurr(ivert) = 1.e-20
          endif
          !
          ! compute absolute differences in wave height and mean period between last 2 iterations
          !
          hsabs = abs ( hscurr(ivert) - hsprev(ivert) )
          tmabs = abs ( tmcurr(ivert) - tmprev(ivert) )
          !
          delhs(ivert) = hsabs
          deltm(ivert) = tmabs
          !
          ! compute curvature of wave height
          !
          hsdif0        = hsdifc(ivert)
          hsdifc(ivert) = 0.5*( hscurr(ivert) - hsprev0 )
          hscurv        = abs ( hsdifc(ivert) - hsdif0 )
          !
          ! compute curvature of mean period
          !
          tmdif0        = tmdifc(ivert)
          tmdifc(ivert) = 0.5*( tmcurr(ivert) - tmprev0 )
          tmcurv        = abs ( tmdifc(ivert) - tmdif0 )
          !
          ! compute required accuracies for wave height
          !
          hsrel  = PNUMS( 1) * hscurr(ivert)
          curvah = PNUMS(15) * hscurr(ivert)
          !
          ! compute required accuracies for mean period
          !
          tmrel  = PNUMS( 1) * tmcurr(ivert)
          curvat = PNUMS(16) * tmcurr(ivert)
          !
          ! count vertices where wave height and period have reached required accuracies
          !
          if ( (hscurv <= curvah .and. hsabs <= max(hsrel,PNUMS(2)) ) .and. &
               (tmcurv <= curvat .and. tmabs <= max(tmrel,PNUMS(3)) ) ) npacct = npacct + 1.
          !
          if (tstfl) then
             if (lhead) write(PRINTF,11)
             write (PRINTF,12) ivert, hsabs, hsabs/hscurr(ivert), hscurv/hscurr(ivert), tmabs, tmabs/tmcurr(ivert), tmcurv/tmcurr(ivert)
             lhead = .false.
          endif
          !
       else
          !
          hscurr(ivert) = 1.e-20
          tmcurr(ivert) = 1.e-20
          !
       endif
       !
    enddo
    !
    ! global sum to npacc and nwetp
    !
    !$omp atomic
    npacc = npacc + npacct
    !$omp atomic
    nwetp = nwetp + nwetpt
    !
!PUN    ! perform global reduction in parallel run
!PUN    !
!PUN    call SwanSumOverNodes ( nwetp )
!PUN    call SwanSumOverNodes ( npacc )
!PUN    !
    ! compute percentage of active vertices where required accuracy has been reached
    !
    !$omp barrier
    !$omp single
    accur = npacc*100./nwetp
    !$omp end single
    !
 11 format(13x,'dHabs          ','dHrel          ','Curvature H    ','dTabs          ','dTrel          ','Curvature T    ')
 12 format(1x,ss,'k=',i7,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
    !
end subroutine SwanConvStopc
