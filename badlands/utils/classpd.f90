!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm class
module classpd

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

    subroutine initialisePD(elevation, demH, pydnodes)

      integer :: pydnodes, k
      real(kind=8),intent(in) :: elevation(pydnodes)
      real(kind=8),intent(inout) :: demH(pydnodes)

      s1 = 0
      s2 = 0

      demH(1:bds) = elevation(1:bds)
      do k = bds+1,pydnodes
          demH(k) = 1.e6
          s1 = s1 + 1
          data1(s1) = k
      enddo

    end subroutine initialisePD

    subroutine fillPD(elevation, sealevel, demH, pydnodes)

      logical :: change
      integer :: pydnodes, n, k, p
      real(kind=8) :: hmin
      real(kind=8),intent(in) :: sealevel
      real(kind=8),intent(in) :: elevation(pydnodes)
      real(kind=8),intent(inout) :: demH(pydnodes)

      change = .true.
      do while(change)
        change = .false.
        do n = 1, s1
          k = data1(n)
          if( demH(k) > elevation(k) )then
            ! Get minimum value
            hmin = 2.e6
            loop: do p = 1, 20
              if( neighbours(k,p) < 0 ) exit loop
              hmin = min(hmin,demH(neighbours(k,p)+1))
            enddo loop
            if( elevation(k) >= hmin + eps )then
              demH(k) = elevation(k)
            else
              if( demH(k) > hmin + eps )then
                demH(k) = hmin + eps
                if(elevation(k)>=sealevel)then
                  if( demH(k) - elevation(k) > fill_TH)then
                    demH(k) = elevation(k) + fill_TH
                  else
                    change = .true.
                  endif
                else
                  if( demH(k) - sealevel > fill_TH)then
                    demH(k) = sealevel + fill_TH
                  else
                    change = .true.
                  endif
                endif
              endif
              s2 = s2 + 1
              data2(s2) = k
            endif
          endif
        enddo
        if( s2 > 0 )then
          s1 = s2
          data1(1:s1) = data2(1:s2)
          s2 = 0
        endif
      enddo

      return

    end subroutine fillPD

    subroutine allfillPD(elevation, demH, pydnodes)

      logical :: change
      integer :: pydnodes, n, k, p
      real(kind=8) :: hmin
      real(kind=8),intent(in) :: elevation(pydnodes)
      real(kind=8),intent(inout) :: demH(pydnodes)

      change = .true.
      do while(change)
        change = .false.
        do n = 1, s1
          k = data1(n)
          if( demH(k) > elevation(k) )then
            ! Get minimum value
            hmin = 2.e6
            loop: do p = 1, 20
              if( neighbours(k,p) < 0 ) exit loop
              hmin = min(hmin,demH(neighbours(k,p)+1))
            enddo loop
            if( elevation(k) >= hmin + eps )then
              demH(k) = elevation(k)
            else
              if( demH(k) > hmin + eps )then
                demH(k) = hmin + eps
                change = .true.
              endif
              s2 = s2 + 1
              data2(s2) = k
            endif
          endif
        enddo
        if( s2 > 0 )then
          s1 = s2
          data1(1:s1) = data2(1:s2)
          s2 = 0
        endif
      enddo

      return

    end subroutine allfillPD

end module classpd
