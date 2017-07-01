subroutine SwanReadfort18
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
!PUN!
!PUN!   Authors
!PUN!
!PUN!   40.95: Marcel Zijlema
!PUN!
!PUN!   Updates
!PUN!
!PUN!   40.95, June 2008: New subroutine
!PUN!
!PUN!   Purpose
!PUN!
!PUN!   Reads fort.18 to obtain the following ADCIRC variables
!PUN!
!PUN!   MNE, MNP, NSTAE, NSTAV, NSTAM and NSTAC
!PUN!
!PUN!   These variables are needed by the message-passing routines in module MESSENGER
!PUN!
!PUN!   Next, we also need the number of elements and nodes in global domain for SWAN computation
!PUN!
!PUN!   Modules used
!PUN!
!PUN    use ocpcomm2
!PUN    use ocpcomm4
!PUN    use SwanGriddata, only: nvertsg, ncellsg
!PUN    use SIZES
!PUN    use GLOBAL, only: NSTAE, NSTAV, NSTAM, NSTAC
!PUN!
!PUN    implicit none
!PUN!
!PUN!   Local variables
!PUN!
!PUN    integer, save           :: ient = 0 ! number of entries in this subroutine
!PUN    integer                 :: idum1    ! dummy integer 1
!PUN    integer                 :: idum2    ! dummy integer 2
!PUN    integer                 :: idum3    ! dummy integer 3
!PUN    integer                 :: idum4    ! dummy integer 4
!PUN    integer                 :: iostat   ! I/O status in call FOR
!PUN    integer                 :: j        ! loop counter
!PUN    character(80)           :: msgfil   ! name of message-passing file including path
!PUN    integer                 :: ndsd     ! unit reference number of file
!PUN    logical                 :: stpnow   ! indicate whether program must be terminated or not
!PUN!
!PUN!   Structure
!PUN!
!PUN!   Description of the pseudo code
!PUN!
!PUN!   Source text
!PUN!
!PUN    if (ltrace) call strace (ient,'SwanReadfort18')
!PUN    !
!PUN    ! open file fort.18
!PUN    !
!PUN    ndsd   = 0
!PUN    iostat = 0
!PUN    msgfil = trim(INPUTDIR)//DIRCH2//'fort.18'
!PUN    call for (ndsd, msgfil, 'OF', iostat)
!PUN    if (stpnow()) goto 900
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, idum3
!PUN    !
!PUN    read(ndsd,'(8x,3i12)', end=950, err=910) ncellsg, idum2, MNE   ! number of elements
!PUN    do j = 1, MNE
!PUN       read(ndsd,'(i12)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) nvertsg, idum2, MNP   ! number of nodes
!PUN    do j = 1, MNP
!PUN       read(ndsd,'(i8)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,'(8x,i8)', end=950, err=910) idum1
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, idum3
!PUN    do j = 1, idum3
!PUN       read(ndsd,'(i8)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAE   ! number of elevation stations
!PUN    do j = 1, NSTAE
!PUN       read(ndsd,'(i8)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAV   ! number of velocity stations
!PUN    do j = 1, NSTAV
!PUN       read(ndsd,'(i8)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAM   ! number of meteorlogical stations
!PUN    do j = 1, NSTAM
!PUN       read(ndsd,'(i8)', end=950, err=910) idum4
!PUN    enddo
!PUN    !
!PUN    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAC   ! number of concentration stations
!PUN    !
!PUN    ! close file fort.18
!PUN    !
!PUN    close(ndsd)
!PUN    !
!PUN 900 return
!PUN    !
!PUN 910 call msgerr (4, 'error reading data from grid file fort.18' )
!PUN    goto 900
!PUN 950 call msgerr (4, 'unexpected end of file in grid file fort.18' )
!PUN    goto 900
!PUN    !
!PUN 100 format((8x,3i8))
!PUN    !
end subroutine SwanReadfort18
