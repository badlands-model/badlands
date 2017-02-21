!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm class
module pdclass

    implicit none

    integer :: dnodes
    integer :: pydx,pydy

    ! Define the data-structure to hold the data
    integer :: s1 = 0
    integer :: s2 = 0
    integer, allocatable, dimension(:) :: data1
    integer, allocatable, dimension(:) :: data2

    ! Set the size of allocated memory blocks
    integer :: block_size

    ! Set neighbourhood array
    integer,allocatable, dimension(:,:) :: neighbours

    ! Set area cells array
    real(kind=8),allocatable, dimension(:) :: area

    integer :: bds, diffnbmax
    real(kind=8) :: eps, fill_TH, diffprop

    ! Sediment distribution parameters
    integer, parameter :: max_it_cyc = 500000
    real(kind=8), parameter :: diff_res = 1.e-2

contains

    subroutine defineparameters

      if(allocated(neighbours)) deallocate(neighbours)
      if(allocated(area)) deallocate(area)
      allocate(neighbours(dnodes,20))
      allocate(area(dnodes))

      if(allocated(data1)) deallocate(data1)
      if(allocated(data2)) deallocate(data2)
      allocate(data1(block_size))
      allocate(data2(block_size))

      return

    end subroutine defineparameters

end module pdclass
