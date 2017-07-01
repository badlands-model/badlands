subroutine SwanReadEasymeshGrid ( basenm, lenfnm )
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
!   Reads Easymesh grid described in <name>.n and <name>.e
!
!   Method
!
!   Grid coordinates of vertices are read from file <name>.n and stored in Swan data structure
!   Vertices of triangles are read from file <name>.e and stored in Swan data structure
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)           :: lenfnm ! length of file names
    character(lenfnm), intent(in) :: basenm ! base name of Easymesh files
!
!   Local variables
!
    character(lenfnm) :: filenm   ! file name
    integer, save     :: ient = 0 ! number of entries in this subroutine
    integer           :: iostat   ! I/O status in call FOR
    integer           :: istat    ! indicate status of allocation
    integer           :: j        ! loop counter
    integer           :: ndsd     ! unit reference number of file
    character(80)     :: line     ! auxiliary textline
    logical           :: stpnow   ! indicate whether program must be terminated or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanReadEasymeshGrid')
    !
    ! open file <name>.n containing the coordinates of vertices
    !
    filenm = trim(basenm)//'.n'
    ndsd   = 0
    iostat = 0
    call for (ndsd, filenm, 'OF', iostat)
    if (stpnow()) goto 900
    !
    ! read first line to determine number of vertices
    !
    read(ndsd, *, end=950, err=910) nverts
    if(.not.allocated(xcugrd)) allocate (xcugrd(nverts), stat = istat)
    if ( istat == 0 ) then
       if(.not.allocated(ycugrd)) allocate (ycugrd(nverts), stat = istat)
    endif
    if ( istat == 0 ) then
       if(.not.allocated(vmark)) allocate (vmark(nverts), stat = istat)
    endif
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanReadEasymeshGrid: array xcugrd, ycugrd or vmark ' )
       goto 900
    endif
    !
    ! read coordinates of vertices and boundary marker
    !
    do j = 1, nverts
       read(ndsd, 100, end=950, err=910) xcugrd(j), ycugrd(j), vmark(j)
    enddo
    !
    ! close file <name>.n
    !
    close(ndsd)
    !
    ! open file <name>.e containing the (Delaunay) triangles
    !
    filenm = trim(basenm)//'.e'
    ndsd   = 0
    iostat = 0
    call for (ndsd, filenm, 'OF', iostat)
    if (stpnow()) goto 900
    !
    ! read first line to determine number of triangles
    !
    read(ndsd, *, end=950, err=910) ncells
    if(.not.allocated(kvertc)) allocate (kvertc(3,ncells), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanReadEasymeshGrid: array kvertc ' )
       goto 900
    endif
    !
    ! read vertices of triangles
    !
    do j = 1, ncells
       read(ndsd, 200, end=950, err=910) kvertc(1,j), kvertc(2,j), kvertc(3,j), line
    enddo
    !
    ! close file <name>.e
    !
    close(ndsd)
    !
    ! Easymesh counters vertices starting from 0 (C style), therefore add 1
    !
    kvertc = kvertc + 1
    !
 900 return
    !
 910 inquire (unit=ndsd, name=filenm)
    call msgerr (4, 'error reading data from Easymesh file '//filenm )
    goto 900
 950 inquire (unit=ndsd, name=filenm)
    call msgerr (4, 'unexpected end of file in Easymesh file '//filenm )
    goto 900
    !
 100 format((6x,2e22.15,i3))
 200 format((5x,3i5,a))
    !
end subroutine SwanReadEasymeshGrid
