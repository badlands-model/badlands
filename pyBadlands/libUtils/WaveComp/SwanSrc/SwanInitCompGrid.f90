subroutine SwanInitCompGrid ( logcom )
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
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   Initialise arrays for description of computational grid
!   in case of unstructured grid
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use m_genarr
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    logical, dimension(6), intent(inout) :: logcom ! give status of which command has been given
!
!   Local variables
!
    integer       :: i        ! loop counter
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: istat    ! indicate status of allocation
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanInitCompGrid')
    !
    ! compute coordinate offsets and reset grid coordinates
    !
    do i = 1, nverts
       if ( .not.LXOFFS ) then
          XOFFS = xcugrd(i)
          YOFFS = ycugrd(i)
          LXOFFS = .true.
          xcugrd(i) = 0.
          ycugrd(i) = 0.
       else
          xcugrd(i) = real(xcugrd(i) - dble(XOFFS))
          ycugrd(i) = real(ycugrd(i) - dble(YOFFS))
       endif
    enddo
    !
    ! check the grid
    !
    call SwanCheckGrid
    !
    ! compute XCGMIN, XCGMAX, YCGMIN, YCGMAX
    !
    XCGMIN =  1.e9
    YCGMIN =  1.e9
    XCGMAX = -1.e9
    YCGMAX = -1.e9
    do i = 1, nverts
       if (xcugrd(i) < XCGMIN) XCGMIN = xcugrd(i)
       if (ycugrd(i) < YCGMIN) YCGMIN = ycugrd(i)
       if (xcugrd(i) > XCGMAX) XCGMAX = xcugrd(i)
       if (ycugrd(i) > YCGMAX) YCGMAX = ycugrd(i)
    enddo
    !
!PUN    XCGMIN = XCGMIN + XOFFS
!PUN    YCGMIN = YCGMIN + YOFFS
!PUN    XCGMAX = XCGMAX + XOFFS
!PUN    YCGMAX = YCGMAX + YOFFS
!PUN    !
!PUN    call SwanMinOverNodes ( XCGMIN )
!PUN    call SwanMinOverNodes ( YCGMIN )
!PUN    call SwanMaxOverNodes ( XCGMAX )
!PUN    call SwanMaxOverNodes ( YCGMAX )
!PUN    !
!PUN    XCGMIN = XCGMIN - XOFFS
!PUN    YCGMIN = YCGMIN - YOFFS
!PUN    XCGMAX = XCGMAX - XOFFS
!PUN    YCGMAX = YCGMAX - YOFFS
!PUN    !
    XCLEN = XCGMAX - XCGMIN
    YCLEN = YCGMAX - YCGMIN
    !
    if(.not.allocated(ac2)) allocate(ac2(MDC,MSC,nverts), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanInitCompGrid: array ac2 ' )
       return
    endif
    ac2 = 0.
    logcom(6) = .true.
    !
    ! set number of vertices and cells in global domain in case of serial run
    !
    if ( nvertsg == 0 ) then
       nvertsg = nverts
       ncellsg = ncells
    endif
    !
    ! the following arrays for structured grids (regular and curvilinear)
    ! are allocated as empty ones
    !
    if ( .not.allocated(KGRPNT) ) allocate(KGRPNT(0,0))
    if ( .not.allocated(KGRBND) ) allocate(KGRBND(0)  )
    !
    ! for sake of convenience, set MCGRD to nverts (for allocating AC1 and COMPDA)
    !
    MCGRD   = nverts
    MCGRDGL = nverts
    !
end subroutine SwanInitCompGrid
