!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements flow parameters computation.
module flowcompute

  implicit none
  
contains

  subroutine discharge(pyStack, pyRcv, pyDischarge, pyDis, pyNodesNb)
  
      integer :: pyNodesNb
      integer,dimension(pyNodesNb),intent(in) :: pyStack
      integer,dimension(pyNodesNb),intent(in) :: pyRcv
      integer,dimension(pyNodesNb),intent(in) :: pyDischarge
      
      integer,dimension(pyNodesNb),intent(out) :: pyDis

      integer :: n, donor, recvr
      
      pyDis = pyDischarge
      
      do n = pyNodesNb, 1, -1
        donor = pyStack(n) + 1
        recvr = pyRcv(donor) + 1
        if( donor /= recvr )then
            pyDis(recvr) = pyDis(recvr) + pyDis(donor)
        endif
      enddo
  
      return

  end subroutine discharge

end module flowcompute
