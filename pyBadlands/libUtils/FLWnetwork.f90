!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module builds the stack used for in Braun & Willett (2013) flow network
! algorithm

module fstack

  use flwstack
  
  implicit none
  
contains

  subroutine build(pyBase,pyRcv,pyDelta,pyDonors,pyStackOrder,pyBaseNb,pyDeltaNb,pyNodesNb)
  
      integer :: pyBaseNb
      integer :: pyDeltaNb
      integer :: pyNodesNb
      integer,dimension(pyBaseNb),intent(in) :: pyBase
      integer,dimension(pyNodesNb),intent(in) :: pyRcv
      integer,dimension(pyDeltaNb),intent(in) :: pyDelta
      
      integer,dimension(pyNodesNb),intent(out) :: pyDonors
      integer,dimension(pyNodesNb),intent(out) :: pyStackOrder

      integer :: p,j,k,success
      integer,dimension(pyNodesNb) :: intArray
      
      j = 0
      
      if(allocated(stackOrder)) deallocate(stackOrder)
      if(allocated(allocs)) deallocate(allocs)
      if(allocated(Donors)) deallocate(Donors)
      if(allocated(Delta)) deallocate(Delta)
      allocate(stackOrder(pyNodesNb))
      allocate(allocs(pyNodesNb))
      allocate(Donors(pyNodesNb))
      allocate(Delta(pyDeltaNb))
      
      intArray = 0
      
      stackOrder = 0
      Delta = pyDelta+1
      
      do k = 1, pyNodesNb
          Donors(Delta(pyRcv(k)+1) + intArray(pyRcv(k)+1)) = k
          intArray(pyRcv(k)+1) = intArray(pyRcv(k)+1)+1
      enddo
      
      allocs = -1
      do p = 1, pyBaseNb
          k = pyBase(p)+1
          j = j+1
          stackOrder(j) = k
          allocs(k) = p
          success = addtostack(p,k,j)
      enddo
  
      pyDonors = Donors-1
      pyStackOrder = stackOrder-1
  
      return

  end subroutine build

end module fstack
