!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm

module pdcompute

  implicit none

contains

  subroutine filling(elevation,pyNgbs,fillTH,epsilon,pybounds,sealimit,demH,pydnodes)

      integer :: pydnodes
      integer,intent(in) :: pybounds
      real(kind=8),intent(in) :: sealimit
      real(kind=8),intent(in) :: fillTH
      real(kind=8),intent(in) :: epsilon
      real(kind=8),intent(in) :: elevation(pydnodes)
      integer,intent(in) :: pyNgbs(pydnodes,20)
      real(kind=8),intent(out) :: demH(pydnodes)

      logical :: flag
      integer :: p, k

      flag = .true.
      demH = 1.e6
      demH(1:pybounds) = elevation(1:pybounds)
      do k=pybounds+1,pydnodes
          if( elevation(k) > sealimit)then
              demH(k) = 1.e6
          else
              demH(k) = elevation(k)
          endif
      enddo

      do while(flag)
        flag=.false.
        do k=1,pydnodes
          if( demH(k) > elevation(k) )then
            loop: do p = 1, 20
                if( pyNgbs(k,p) < 0) exit loop
                if( elevation(k) >= demH(pyNgbs(k,p)+1) + epsilon )then
                    demH(k) = elevation(k)
                else
                    if( demH(k) > demH(pyNgbs(k,p)+1) + epsilon )then
                        demH(k) = demH(pyNgbs(k,p)+1) + epsilon
                        if(demH(k) - elevation(k) > fillTH)then
                            demH(k) = elevation(k) + fillTH
                        else
                            flag=.true.
                        endif
                    endif
                endif
            enddo loop
          endif
        enddo
      enddo

      return

  end subroutine filling

end module pdcompute

module pdstack

  implicit none

  ! Define the data-structure to hold the data
  integer :: s1 = 0
  integer :: s2 = 0
  integer, allocatable, dimension(:) :: data1
  integer, allocatable, dimension(:) :: data2

  ! Set the size of allocated memory blocks
  integer :: block_size

  ! Set neighbourhood array
  integer,allocatable, dimension(:,:) :: neighbours

  integer :: bds
  real(kind=8) :: eps, fill_TH

contains

  subroutine initialisePD(elevation, demH, sealimit, pydnodes)

    integer :: pydnodes, k
    real(kind=8),intent(in) :: sealimit
    real(kind=8),intent(in) :: elevation(pydnodes)
    real(kind=8),intent(inout) :: demH(pydnodes)

    s1 = 0
    s2 = 0

    demH(1:bds) = elevation(1:bds)
    do k = bds+1,pydnodes
        if( elevation(k) > sealimit)then
            demH(k) = 1.e6
            s1 = s1 + 1
            data1(s1) = k
        else
            demH(k) = elevation(k)
        endif
    enddo

  end subroutine initialisePD

  subroutine fillPD(elevation, demH, pydnodes)

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
              if( demH(k) - elevation(k) > fill_TH )then
                demH(k) = elevation(k) + fill_TH
              else
                change = .true.
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

  subroutine pitparams(pyNgbs,fillTH,epsilon,pybounds,pydnodes)

    integer :: pydnodes
    integer,intent(in) :: pybounds
    real(kind=8),intent(in) :: fillTH
    real(kind=8),intent(in) :: epsilon
    integer,intent(in) :: pyNgbs(pydnodes,20)

    if(allocated(neighbours)) deallocate(neighbours)
    allocate(neighbours(pydnodes,20))
    neighbours = pyNgbs
    bds = pybounds
    block_size = pydnodes - bds
    eps = epsilon
    fill_TH = fillTH
    if(allocated(data1)) deallocate(data1)
    if(allocated(data2)) deallocate(data2)
    allocate(data1(block_size))
    allocate(data2(block_size))

    return

  end subroutine pitparams

  subroutine pitfilling(elevation,sealimit,allfill,demH,pydnodes)

    integer :: pydnodes
    integer,intent(in) :: allfill
    real(kind=8),intent(in) :: sealimit
    real(kind=8),intent(in) :: elevation(pydnodes)
    real(kind=8),intent(out) :: demH(pydnodes)

    ! Initialisation phase
    call initialisePD(elevation,demH,sealimit,pydnodes)

    ! Filling phase
    if(allfill == 0)then
      call fillPD(elevation,demH,pydnodes)
    else
      call allfillPD(elevation,demH,pydnodes)
    endif

    return

  end subroutine pitfilling

end module pdstack
