subroutine SwanThreadBounds ( nwetp, ivlow, ivup, tlist )
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
!   41.10: Marcel Zijlema
!
!   Updates
!
!   41.10, August 2009: New subroutine
!
!   Purpose
!
!   Determines load-balanced loop bounds for calling thread
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(out)                    :: ivlow     ! lower index in range of vertices in calling thread
    integer, intent(out)                    :: ivup      ! upper index in range of vertices in calling thread
    integer, dimension(nverts), intent(out) :: tlist     ! vertex list for calling thread
    !
    real, intent(in)                        :: nwetp     ! total number of active vertices
!
!   Local variables
!
    integer                                 :: i         ! counter
    integer, save                           :: ient = 0  ! number of entries in this subroutine
    integer                                 :: ith       ! thread counter
    integer                                 :: ivert     ! vertex index
    integer                                 :: j         ! counter
    integer                                 :: kvert     ! loop counter over vertices
    integer, dimension(10)                  :: nacvt     ! number of active vertices for i-th thread
    integer                                 :: ncurvt    ! number of currently assigned vertices to a thread
    integer                                 :: nvcum     ! cumulative number of vertices
    integer                                 :: nth       ! number of threads
    integer                                 :: tid       ! thread number
    !
    type(verttype), dimension(:), pointer   :: vert      ! datastructure for vertices with their attributes
    !
!$  integer, external                       :: omp_get_num_threads ! number of OpenMP threads being used
!$  integer, external                       :: omp_get_thread_num  ! get thread number
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanThreadBounds')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    ! determine number of threads and thread number
    !
    nth = 1
    tid = 0
    !$ nth = omp_get_num_threads()
    !$ tid = omp_get_thread_num()
    tid = tid + 1
    !
    ! determine load-balanced sizes for all threads
    !
    nvcum = 0
    do i = 1, nth
       nacvt(i) = (nint(nwetp)*i)/nth - nvcum
       nvcum    = (nint(nwetp)*i)/nth
    enddo
    !
    ! determine loop bounds for calling thread
    !
    ivlow = nverts+1
    ivup  = 0
    !
    ith    = 1
    ncurvt = 0
    !
    do kvert = 1, nverts
       !
       ivert = vlist(kvert)
       !
       if ( vert(ivert)%active ) then
          !
          if ( ith == tid ) then
             ivlow = min(kvert,ivlow)
             ivup  = max(kvert,ivup )
          endif
          ncurvt = ncurvt + 1
          !
          if ( ncurvt >= nacvt(ith) ) then
             ith    = ith + 1
             ncurvt = 0
          endif
          !
       endif
       !
    enddo
    !
    ! determine vertex list for calling thread based on vlist
    ! first active vertices followed by inactive ones
    !
    i = ivlow
    j = ivup
    !
    do kvert = ivlow, ivup
       !
       ivert = vlist(kvert)
       !
       if ( vert(ivert)%active ) then
          !
          tlist(i) = ivert
          i = i + 1
          !
       else
          !
          tlist(j) = ivert
          j = j - 1
          !
       endif
       !
    enddo
    !
end subroutine SwanThreadBounds
