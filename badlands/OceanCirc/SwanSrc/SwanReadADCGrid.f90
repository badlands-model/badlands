subroutine SwanReadADCGrid
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
!   40.95: Marcel Zijlema
!   41.07: Casey Dietrich
!
!   Updates
!
!   40.80, December 2007: New subroutine
!   40.95,     June 2008: parallelization of unSWAN using MESSENGER of ADCIRC
!   41.07,   August 2009: use ADCIRC boundary info to mark all boundary vertices
!
!   Purpose
!
!   Reads ADCIRC grid described in fort.14
!
!   Method
!
!   Grid coordinates of vertices are read from file fort.14 and stored in Swan data structure
!   Vertices of triangles are read from file fort.14 and stored in Swan data structure
!
!   Bottom topography from file fort.14 will also be stored
!
!   Modules used
!
!PUN    use ocpcomm2
    use ocpcomm4
    use m_genarr
    use SwanGriddata
!PUN    use SIZES
!PUN    use MESSENGER
!
    implicit none
!
!   Local variables
!
    character(80)           :: grdfil   ! name of grid file including path
    integer, save           :: ient = 0 ! number of entries in this subroutine
    integer                 :: idum     ! dummy integer
    integer                 :: ii       ! auxiliary integer
    integer                 :: iostat   ! I/O status in call FOR
    integer                 :: istat    ! indicate status of allocation
    integer                 :: itype    ! ADCIRC boundary type
    integer                 :: ivert    ! vertex index
    integer                 :: ivert1   ! another vertex index
    integer                 :: j        ! loop counter
    integer                 :: k        ! loop counter
    integer                 :: n1       ! auxiliary integer
    integer                 :: n2       ! another auxiliary integer
    integer                 :: ndsd     ! unit reference number of file
    integer                 :: nopbc    ! number of open boundaries in ADCIRC
    integer                 :: vm       ! boundary marker
    character(80)           :: line     ! auxiliary textline
    logical                 :: stpnow   ! indicate whether program must be terminated or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanReadADCGrid')
    !
    ! open file fort.14
    !
    ndsd   = 0
    iostat = 0
    grdfil = 'fort.14'
!PUN    grdfil = trim(INPUTDIR)//DIRCH2//trim(grdfil)
    call for (ndsd, grdfil, 'OF', iostat)
    if (stpnow()) goto 900
    !
    ! skip first line
    !
    read(ndsd,'(a80)', end=950, err=910) line
    !
    ! read number of elements and number of vertices
    !
    read(ndsd, *, end=950, err=910) ncells, nverts
    if(.not.allocated(xcugrd)) allocate (xcugrd(nverts), stat = istat)
    if ( istat == 0 ) then
       if(.not.allocated(ycugrd)) allocate (ycugrd(nverts), stat = istat)
    endif
    if ( istat == 0 ) then
       if(.not.allocated(DEPTH)) allocate (DEPTH(nverts), stat = istat)
    endif
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanReadADCGrid: array xcugrd, ycugrd or depth ' )
       goto 900
    endif
    !
    ! read coordinates of vertices and bottom topography
    !
    do j = 1, nverts
       read(ndsd, *, end=950, err=910) ii, xcugrd(ii), ycugrd(ii), DEPTH(ii)
       if ( ii/=j ) call msgerr ( 1, 'numbering of vertices is not sequential in grid file fort.14 ' )
    enddo
    !
    if(.not.allocated(kvertc)) allocate (kvertc(3,ncells), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanReadADCGrid: array kvertc ' )
       goto 900
    endif
    !
    ! read vertices of triangles
    !
    do j = 1, ncells
       read(ndsd, *, end=950, err=910) ii, idum, kvertc(1,ii), kvertc(2,ii), kvertc(3,ii)
       if ( ii/=j ) call msgerr ( 1, 'numbering of triangles is not sequential in grid file fort.14 ' )
    enddo
    !
    if(.not.allocated(vmark)) allocate (vmark(nverts), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanReadADCGrid: array vmark ' )
       goto 900
    endif
    vmark = 0
    !
    ! read ADCIRC boundary information and store boundary markers
    !
    read(ndsd, *, end=950, err=910) nopbc
    read(ndsd, *, end=950, err=910) idum
    do j = 1, nopbc
!PUN       if ( MNPROC==1 ) then
          vm = j
          read(ndsd, *, end=950, err=910) n2
!PUN       else
!PUN          read(ndsd, *, end=950, err=910) n2, vm
!PUN       endif
       do k = 1, n2
           read(ndsd, *, end=950, err=910) ivert
           vmark(ivert) = vm
       enddo
    enddo
    !
    read(ndsd, *, end=950, err=910) n1
    read(ndsd, *, end=950, err=910) idum
    do j = 1, n1
!PUN       if ( MNPROC==1 ) then
          vm = nopbc + j
          read(ndsd, *, end=950, err=910) n2, itype
!PUN       else
!PUN          read(ndsd, *, end=950, err=910) n2, itype, vm
!PUN       endif
       if ( itype /= 4 .and. itype /= 24 ) then
          do k = 1, n2
             read(ndsd, *, end=950, err=910) ivert
             vmark(ivert) = vm
          enddo
       else
          do k = 1, n2
             read(ndsd, *, end=950, err=910) ivert, ivert1
             vmark(ivert ) = vm
             vmark(ivert1) = vm
          enddo
       endif
    enddo
    !
    ! close file fort.14
    !
    close(ndsd)
    !
!PUN       ! ghost vertices are marked with +999
!PUN       !
!PUN       do j = 1, NEIGHPROC
!PUN          do k = 1, NNODRECV(j)
!PUN             ivert = IRECVLOC(k,j)
!PUN             vmark(ivert) = 999
!PUN          enddo
!PUN       enddo
!PUN       !
 900 return
    !
 910 call msgerr (4, 'error reading data from grid file fort.14' )
    goto 900
 950 call msgerr (4, 'unexpected end of file in grid file fort.14' )
    goto 900
    !
end subroutine SwanReadADCGrid
