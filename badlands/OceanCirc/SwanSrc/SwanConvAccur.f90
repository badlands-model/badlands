subroutine SwanConvAccur ( accur, hscurr, tmcurr, delhs, deltm, xytst, spcsig, ac2, ivlow, ivup )
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
!   41.10: Marcel Zijlema
!
!   Updates
!
!   40.80, October 2007: New subroutine
!   41.10,  August 2009: parallelization using OpenMP directives
!
!   Purpose
!
!   Determine accuracy of some integral parameters for convergence check
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
    real, dimension(nverts), intent(inout)      :: tmcurr ! mean period at current iteration level
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
    real                                  :: fact     ! auxiliary factor
    real                                  :: hsabs    ! absolute difference in wave height between last 2 iterations
    real                                  :: hsoval   ! required accuracy with respect to overall error in wave height
    real                                  :: hsprev   ! wave height at previous iteration level
    real                                  :: hsrel    ! required accuracy with respect to relative error in wave height
    real                                  :: hsmean   ! space-averaged wave height at previous iteration level
    real                                  :: hsmeant  ! hsmean variable for calling thread
    real                                  :: m0       ! moment of zeroth order
    real                                  :: m1       ! moment of first order
    real                                  :: npacc    ! number of vertices in which required accuracy has been reached
    real                                  :: npacct   ! npacc counter for calling thread
    real                                  :: nwetp    ! total number of active vertices
    real                                  :: nwetpt   ! nwetp counter for calling thread
    real                                  :: tmmean   ! space-averaged mean period at previous iteration level
    real                                  :: tmmeant  ! tmmean variable for calling thread
    real                                  :: tmabs    ! absolute difference in mean period between last 2 iterations
    real                                  :: tmoval   ! required accuracy with respect to overall error in mean period
    real                                  :: tmprev   ! mean period at previous iteration level
    real                                  :: tmrel    ! required accuracy with respect to relative error in mean period
    !
    logical                               :: lhead    ! logical indicating to write header
    logical                               :: tstfl    ! indicates whether vertex is a test point
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!   Common variables
!
    common/convaccur/npacc,nwetp,hsmean,tmmean
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanConvAccur')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    !$omp single
    npacc = 0.
    nwetp = 0.
    !
    hsmean = 0.
    tmmean = 0.
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
    ! calculate space-averaged wave height and mean period
    !
    hsmeant = 0.
    tmmeant = 0.
    !
    do ivert = ivlow, ivup
       !
       if ( vert(ivert)%active ) then
          !
          nwetpt = nwetpt + 1.
          !
          hsmeant = hsmeant + hscurr(ivert)
          tmmeant = tmmeant + tmcurr(ivert)
          !
       endif
       !
    enddo
    !
    ! global sum to nwetp, hsmean and tmmean
    !
    !$omp atomic
    nwetp  = nwetp + nwetpt
    !$omp atomic
    hsmean = hsmean + hsmeant
    !$omp atomic
    tmmean = tmmean + tmmeant
    !
!PUN    ! perform global reductions in parallel run
!PUN    !
!PUN    call SwanSumOverNodes ( nwetp  )
!PUN    call SwanSumOverNodes ( hsmean )
!PUN    call SwanSumOverNodes ( tmmean )
!PUN    !
    !$omp barrier
    !$omp single
    hsmean = hsmean/nwetp
    tmmean = tmmean/nwetp
    !$omp end single
    !
    ! calculate a set of accuracy parameters based on relative and overall accuracy measures for Hs and Tm
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
          ! store wave height and mean period of previous iteration level
          !
          hsprev = max( 1.e-20, hscurr(ivert) )
          tmprev = max( 1.e-20, tmcurr(ivert) )
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
          hsabs = abs ( hscurr(ivert) - hsprev )
          tmabs = abs ( tmcurr(ivert) - tmprev )
          !
          delhs(ivert) = hsabs
          deltm(ivert) = tmabs
          !
          ! compute required accuracies for wave height and mean period
          !
          hsrel  = PNUMS( 1) * hsprev
          hsoval = PNUMS(15) * hsmean
          !
          tmrel  = PNUMS( 1) * tmprev
          tmoval = PNUMS(16) * tmmean
          !
          ! count vertices where wave height and mean period have reached required accuracies
          !
          if ( hsabs <= max(hsrel,hsoval) .and. tmabs <= max(tmrel,tmoval) ) npacct = npacct + 1.
          !
          if (tstfl) then
             if (lhead) write(PRINTF,11)
             write (PRINTF,12) ivert, hsabs/hsprev, hsabs/hsmean, tmabs/tmprev, tmabs/tmmean
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
    ! global sum to npacc
    !
    !$omp atomic
    npacc = npacc + npacct
    !
!PUN    ! perform global reduction in parallel run
!PUN    !
!PUN    call SwanSumOverNodes ( npacc )
!PUN    !
    ! compute percentage of active vertices where required accuracy has been reached
    !
    !$omp barrier
    !$omp single
    accur = npacc*100./nwetp
    !$omp end single
    !
 11 format(13x,'dHrel          ','dHoval         ','dTm01rel       ','dTm01oval ')
 12 format(1x,ss,'k=',i7,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
    !
end subroutine SwanConvAccur
