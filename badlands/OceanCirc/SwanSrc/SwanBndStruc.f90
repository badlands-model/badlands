subroutine SwanBndStruc ( xcgrid, ycgrid )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.14: Nico Booij
!
!   Updates
!
!   41.14, Aug. 2010: New subroutine
!
!   Purpose
!
!   Generates output curves 'BOUNDARY' and 'BOUND_01' until 'BOUND_04'
!   in the case of a structured grid (both regular and curvilinear).
!
!   Method
!
!   The procedure starts at the origin of the grid, and then moves
!   around the grid in counterclockwise order.
!
!   Modules used
!
    use ocpcomm4
    use SWCOMM1
    use SWCOMM2
    use SWCOMM3
    use OUTP_DATA
!
    implicit none
!
!   Argument variables
!
    real, dimension(MXC,MYC), intent(in) :: xcgrid   ! x-coordinate of computational grid
    real, dimension(MXC,MYC), intent(in) :: ycgrid   ! y-coordinate of computational grid
!
!   Local variables
!
    integer                              :: ix, iy   ! point index
    integer                              :: ibnd     ! number of valid points along boundary
    integer                              :: iside    ! side counter (1..4)
    integer                              :: lside    ! number of points on one side
    integer                              :: ispt     ! number of valid points on one side
    integer                              :: mip      ! number of output points
    integer                              :: xstep, ystep
    integer                              :: ix0, ix1, iy0, iy1
    integer                              :: ii, jj   ! counters
    integer, save                        :: ient = 0 ! number of entries in this subroutine
    !
    logical, save                        :: done=.false. ! if true procedure has been done
    logical                              :: EQREAL   ! function
    !
    real                                 :: xp, yp       ! one boundary point
    real, allocatable, dimension (:)     :: xbnd, ybnd   ! points of whole boundary
    real, allocatable, dimension (:)     :: xsid, ysid   ! points of one side
    !
    character(80)                        :: msgstr   ! string to pass message
    character(len=8)                     :: psname   ! name assigned to output curve
    !
    TYPE(OPSDAT), POINTER :: OPSTMP, ROPS
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanBndStruc')
    !
    ! if list of boundary vertices is already filled, return
    !
    if (done) return
    !
    allocate (xbnd(1:2*(mxc+myc-2)))
    allocate (ybnd(1:2*(mxc+myc-2)))
    ibnd = 0
    do iside = 1, 4
      if (iside==1) then
         ix0 = 1
         ix1 = mxc
         xstep = 1
         iy0 = 1
         iy1 = 1
         ystep = 0
         lside = mxc
      elseif (iside==2) then
         ix0 = mxc
         ix1 = mxc
         xstep = 0
         iy0 = 1
         iy1 = myc
         ystep = 1
         lside = myc
      elseif (iside==3) then
         ix0 = mxc
         ix1 = 1
         xstep = -1
         iy0 = myc
         iy1 = myc
         ystep = 0
         lside = mxc
      elseif (iside==4) then
         ix0 = 1
         ix1 = 1
         xstep = 0
         iy0 = myc
         iy1 = 1
         ystep = -1
         lside = myc
      endif
      !
      allocate (xsid(1:lside))
      allocate (ysid(1:lside))
      ix = ix0
      iy = iy0
      ispt = 0
      !
      do ii = 1, lside
        !
        ! loop over points of one side of the grid
        !
        xp = xcgrid(ix,iy)
        yp = ycgrid(ix,iy)
        if (.not.(EQREAL(xp,OVEXCV(1)).or.EQREAL(yp,OVEXCV(2)))) then
          ! point has valid coordinates
          ispt = ispt + 1
          xsid(ispt) = xp
          ysid(ispt) = yp
          if (ii<lside) then
            ! corner point appears only once in curve BOUNDARY
            ibnd = ibnd + 1
            xbnd(ibnd) = xp
            ybnd(ibnd) = yp
          endif
        endif
        if (ystep==0) then
          ix = ix + xstep
        else
          iy = iy + ystep
        endif
      enddo
      !
      write (psname, '(A6, I2.2)') 'BOUND_', iside
      mip = ispt
      if (mip>0) then
        allocate(OPSTMP)
        OPSTMP%PSNAME = psname
        OPSTMP%PSTYPE = 'C'
        OPSTMP%MIP = mip
        allocate(OPSTMP%XP(mip))
        allocate(OPSTMP%YP(mip))
        do jj = 1, mip
           OPSTMP%XP(jj) = xsid(jj)
           OPSTMP%YP(jj) = ysid(jj)
        enddo
        deallocate (xsid, ysid)
        nullify (OPSTMP%NEXTOPS)
        if ( .not.LOPS ) then
           FOPS = OPSTMP
           COPS => FOPS
           LOPS = .TRUE.
        else
           COPS%NEXTOPS => OPSTMP
           COPS => OPSTMP
        endif
        if (ITEST>=10) write (PRTEST, *) 'Output curve ', psname,    &
                       ' with ', mip, ' points is generated'
      else
        call MSGERR(1,'No output points found in '//psname)
      endif
    enddo
    !
    psname = 'BOUNDARY'
    mip = ibnd
    if (mip>0) then
      allocate(OPSTMP)
      OPSTMP%PSNAME = psname
      OPSTMP%PSTYPE = 'C'
      OPSTMP%MIP = mip
      allocate(OPSTMP%XP(mip))
      allocate(OPSTMP%YP(mip))
      do jj = 1, mip
         OPSTMP%XP(jj) = xbnd(jj)
         OPSTMP%YP(jj) = ybnd(jj)
      enddo
      deallocate(xbnd,ybnd)
      nullify (OPSTMP%NEXTOPS)
      if ( .not.LOPS ) then
         FOPS = OPSTMP
         COPS => FOPS
         LOPS = .TRUE.
      else
         COPS%NEXTOPS => OPSTMP
         COPS => OPSTMP
      endif
      if (ITEST>=10) write (PRTEST, *) 'Output curve ', psname,    &
                     ' with ', mip, ' points is generated'
    else
      call MSGERR(1,'No output points found in '//psname)
    endif
    !
    done = .true.  ! prevents second entry into this subroutine
    !
end subroutine SwanBndStruc
