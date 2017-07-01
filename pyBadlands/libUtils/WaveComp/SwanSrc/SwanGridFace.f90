subroutine SwanGridFace ( nfaces, ncells, nverts, xcugrd, ycugrd, kvertf )
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
!   Fills face-based data structure
!
!   Method
!
!   Based on unstructured grid
!
!   Modules used
!
    use ocpcomm4
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                       :: ncells  ! number of cells in grid
    integer, intent(in)                       :: nfaces  ! number of faces in grid
    integer, intent(in)                       :: nverts  ! number of vertices in grid
    !
    integer, dimension(2, nfaces), intent(in) :: kvertf  ! vertices of the face
                                                         ! (must be filled by a gridgenerator!)
    !
    real, dimension(nverts), intent(in)       :: xcugrd  ! the x-coordinates of the grid vertices
    real, dimension(nverts), intent(in)       :: ycugrd  ! the y-coordinates of the grid vertices
!
!   Local variables
!
    integer                               :: icell     ! loop counter over cells
    integer                               :: icell1    ! sequence number of cell 1 adjacent to present face
    integer                               :: icell2    ! sequence number of cell 2 adjacent to present face
    integer                               :: iface     ! loop counter over faces
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: ivert     ! loop counter over vertices
    integer                               :: j         ! loop counter
    integer                               :: jf        ! loop counter
    integer                               :: k         ! counter
    integer                               :: v1        ! first vertex of present face
    integer                               :: vl1       ! first vertex of local face for given cell
    integer                               :: v2        ! second vertex of present face
    integer                               :: vl2       ! second vertex of local face for given cell
    !
    integer, dimension(:  ), allocatable  :: cntv1     ! array of vertex counter for for vertex 1
    integer, dimension(:  ), allocatable  :: cntv2     ! array of vertex counter for for vertex 2
    integer, dimension(:,:), allocatable  :: iflist1   ! list of index faces stored for vertex 1
    integer, dimension(:,:), allocatable  :: iflist2   ! list of index faces stored for vertex 2
    !
    real                                  :: carea1    ! area of cell 1 adjacent to present face
    real                                  :: carea2    ! area of cell 2 adjacent to present face
    real                                  :: lengthf   ! length of present face
    real                                  :: xdiff     ! difference in x-coordinate between vertex 2 and vertex 1
    real                                  :: ydiff     ! difference in y-coordinate between vertex 2 and vertex 1
    !
    logical                               :: facefound ! true if face is found
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face      ! datastructure for faces with their attributes
    type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGridFace')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! loop over all faces
    !
    do iface = 1, nfaces
       !
       ! identification number
       !
       face(iface)%atti(FACEID) = iface
       !
       ! Fill vertices
       !
       v1 = kvertf(1,iface)
       v2 = kvertf(2,iface)
       face(iface)%atti(FACEV1) = v1
       face(iface)%atti(FACEV2) = v2
       !
       ! compute length of face
       !
       xdiff   = xcugrd(v2) - xcugrd(v1)
       ydiff   = ycugrd(v2) - ycugrd(v1)
       lengthf = sqrt( xdiff*xdiff + ydiff*ydiff )
       !
       ! fill length of face and normal to face
       !
       face(iface)%attr(FACELEN  ) = lengthf
       face(iface)%attr(FACENORMX) =  ydiff/lengthf
       face(iface)%attr(FACENORMY) = -xdiff/lengthf
       !
       ! compute coordinates of midpoint of face
       !
       face(iface)%attr(FACEMX) = 0.5*(xcugrd(v1) + xcugrd(v2))
       face(iface)%attr(FACEMY) = 0.5*(ycugrd(v1) + ycugrd(v2))
       !
       face(iface)%atti(FACEC1 ) = 0
       face(iface)%atti(FACEC2 ) = 0
       face(iface)%atti(FMARKER) = 0
       !
    enddo
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
    do iface = 1, nfaces
       !
       v1 = face(iface)%atti(FACEV1)
       v2 = face(iface)%atti(FACEV2)
       !
       k = cntv1(v1) +1
       if ( k > 10 ) then
          call msgerr ( 4, 'SwanGridFace: more than 10 faces around vertex ' )
          return
       endif
       cntv1  (v1  ) = k
       iflist1(v1,k) = iface
       !
       k = cntv2(v2) +1
       if ( k > 10 ) then
          call msgerr ( 4, 'SwanGridFace: more than 10 faces around vertex ' )
          return
       endif
       cntv2  (v2  ) = k
       iflist2(v2,k) = iface
       !
    enddo
    !
    ! loop over all cells
    !
    do icell = 1, ncells
       !
       ! loop over all local faces of the cell
       !
       do jf = 1, cell(icell)%nof
          !
          ! determine vertices of the local face
          !
          vl1 = cell(icell)%face(jf)%atti(FACEV1)
          vl2 = cell(icell)%face(jf)%atti(FACEV2)
          !
          ! search for identification number of that face
          !
          facefound = .false.
          !
          kloop: do k = 1, 10
             !
             iface = iflist1(vl1,k)
             !
             do j = 1, 10
                if ( iflist2(vl2,j) == iface ) then
                   facefound = .true.
                   exit kloop
                endif
             enddo
             !
          enddo kloop
          !
          if ( .not.facefound ) then
             !
             jloop: do j = 1, 10
                !
                iface = iflist2(vl1,j)
                !
                do k = 1, 10
                   if ( iflist1(vl2,k) == iface ) then
                      facefound = .true.
                      exit jloop
                   endif
                enddo
                !
             enddo jloop
             !
          endif
          !
          if ( facefound ) then
             cell(icell)%face(jf)%atti(FACEID) = iface
          else
             call msgerr ( 4, 'inconsistency found in SwanGridFace: no face found ' )
             return
          endif
          !
          v1 = face(iface)%atti(FACEV1)
          v2 = face(iface)%atti(FACEV2)
          !
          ! Requirement: face(iface)%atti(FACEC1) < face(iface)%atti(FACEC2)
          !
          if ( v1 == vl1 ) then
             if ( face(iface)%atti(FACEC1) == 0 ) then
                face(iface)%atti(FACEC1) = icell
             else
                call msgerr ( 4, 'SwanGridFace: not all cells have counterclockwise order of vertices ' )
                return
             endif
          else
             if ( face(iface)%atti(FACEC2) == 0 ) then
                face(iface)%atti(FACEC2) = icell
             else
                call msgerr ( 4, 'SwanGridFace: not all cells have counterclockwise order of vertices ' )
                return
             endif
          endif
          !
       enddo
       !
    enddo
    !
    deallocate(cntv1,cntv2,iflist1,iflist2)
    !
    ! loop over all faces
    !
    do iface = 1, nfaces
       !
       ! marks boundary face
       !
       if (face(iface)%atti(FACEC2) == 0) face(iface)%atti(FMARKER) = 1
       !
       ! store interpolation factors
       !
       if ( face(iface)%atti(FMARKER) == 1 ) then
          !
          face(iface)%attr(FACELINPF) = 0.
          !
       else
          !
          icell1 = face(iface)%atti(FACEC1)
          icell2 = face(iface)%atti(FACEC2)
          !
          carea1 = cell(icell1)%attr(CELLAREA)
          carea2 = cell(icell2)%attr(CELLAREA)
          !
          face(iface)%attr(FACELINPF) = carea1/(carea1+carea2)
          !
       endif
       !
    enddo
    !
    ! marks boundary cells
    !
    do icell = 1, ncells
       !
       cell(icell)%atti(CMARKER) = 0
       !
       ! loop over all local faces of the cell
       !
       do jf = 1, cell(icell)%nof
          !
          iface = cell(icell)%face(jf)%atti(FACEID)
          !
          if ( face(iface)%atti(FMARKER) == 1 ) then
             cell(icell)%atti(CMARKER) = 1
             exit
          endif
          !
       enddo
       !
    enddo

end subroutine SwanGridFace
