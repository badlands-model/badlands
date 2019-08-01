subroutine SwanCreateEdges
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
!   Generates edge-based data structure
!
!   Method
!
!   Faces of triangles are computed from elements and stored in Swan data structure
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
!
    implicit none
!
!   Local variables
!
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: iface     ! actual face of present cell
    integer                               :: istat     ! indicate status of allocation
    integer                               :: j         ! loop counter
    integer                               :: k         ! loop counter
    integer                               :: m         ! counter
    integer                               :: maxnf     ! maximum number of faces
    integer                               :: n         ! counter
    integer                               :: v1        ! first vertex of present cell
    integer                               :: v2        ! second vertex of present cell
    integer                               :: v3        ! third vertex of present cell
    integer, dimension(2,3)               :: vf        ! vertices of faces of present cell
    !
    integer, dimension(:  ), allocatable  :: cntv1     ! array of vertex counter for for vertex 1
    integer, dimension(:  ), allocatable  :: cntv2     ! array of vertex counter for for vertex 2
    integer, dimension(:,:), allocatable  :: iflist1   ! list of index faces stored for vertex 1
    integer, dimension(:,:), allocatable  :: iflist2   ! list of index faces stored for vertex 2
    !
    logical                               :: facefound ! true if face is found
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanCreateEdges')
    !
    ! determine and store vertices of faces
    !
    maxnf = 3*ncells   ! maximum number of faces
    !
    if(.not.allocated(kvertf)) allocate (kvertf(2,maxnf), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanCreateEdges: array kvertf ' )
       return
    endif
    kvertf = 0
    !
    allocate(cntv1  (nverts   ))
    allocate(cntv2  (nverts   ))
    allocate(iflist1(nverts,10))
    allocate(iflist2(nverts,10))
    !
    cntv1   =  0
    cntv2   =  0
    iflist1 = -1
    iflist2 = -2
    !
    nfaces = 0
    !
    do j = 1, ncells
       !
       v1 = kvertc(1,j)
       v2 = kvertc(2,j)
       v3 = kvertc(3,j)
       !
       vf(1,1) = v2
       vf(2,1) = v3
       vf(1,2) = v3
       vf(2,2) = v1
       vf(1,3) = v1
       vf(2,3) = v2
       !
       kloop: do k = 1, 3
          !
          ! get two vertices of a face
          !
          v1 = vf(1,k)
          v2 = vf(2,k)
          !
          ! search for identification number of that face
          !
          facefound = .false.
          !
          mloop: do m = 1, 10
             !
             iface = iflist1(v1,m)
             !
             do n = 1, 10
                if ( iflist2(v2,n) == iface ) then
                   facefound = .true.
                   exit mloop
                endif
             enddo
             !
          enddo mloop
          !
          if ( .not.facefound ) then
             !
             nloop: do n = 1, 10
                !
                iface = iflist2(v1,n)
                !
                do m = 1, 10
                   if ( iflist1(v2,m) == iface ) then
                      facefound = .true.
                      exit nloop
                   endif
                enddo
                !
             enddo nloop
             !
          endif
          !
          if ( facefound ) cycle kloop
          !
          nfaces = nfaces + 1
          if ( nfaces > maxnf ) then
             call msgerr ( 4, 'inconsistency found in SwanCreateEdges: more than maximum allowable faces found ' )
             return
          endif
          !
          m = cntv1(v1) +1
          if ( m > 10 ) then
             call msgerr ( 4, 'SwanCreateEdges: more than 10 faces around vertex ' )
             return
          endif
          cntv1  (v1  ) = m
          iflist1(v1,m) = nfaces
          !
          m = cntv2(v2) +1
          if ( m > 10 ) then
             call msgerr ( 4, 'SwanCreateEdges: more than 10 faces around vertex ' )
             return
          endif
          cntv2  (v2  ) = m
          iflist2(v2,m) = nfaces
          !
          kvertf(1,nfaces) = v1
          kvertf(2,nfaces) = v2
          !
       enddo kloop
       !
    enddo
    !
    deallocate(cntv1,cntv2,iflist1,iflist2)
    !
end subroutine SwanCreateEdges
