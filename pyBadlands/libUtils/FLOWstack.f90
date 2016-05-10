!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements the stack used for in Braun & Willett (2013) flow network
! algorithm

module flwstack

  implicit none

  integer,dimension(:),allocatable :: allocs
  integer,dimension(:),allocatable :: Donors
  integer,dimension(:),allocatable :: Delta
  integer,dimension(:),allocatable :: stackOrder
  
contains

  recursive function addtostack(base,donor,stackID) result(success)
  
      integer :: base,donor,stackID,n,success

      success = 1
      do n = Delta(donor),Delta(donor+1)-1
          if((allocs(Donors(n)) /= base))then
              donor = Donors(n)
              stackID = stackID+1
              stackOrder(stackID) = donor
              allocs(donor) = base
              success = addtostack(base,donor,stackID)
          endif
      enddo
      success = 0

  end function addtostack

end module flwstack