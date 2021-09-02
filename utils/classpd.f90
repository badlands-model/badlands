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
    real(kind=8) :: eps, fill_TH, diffprop, propA, propB

    ! Sediment distribution parameters
    integer, parameter :: max_it_cyc = 500000
    real(kind=8), parameter :: diff_res = 1.e-2

    ! Queue node definition: index and elevation
    type node
      integer :: id
      real(kind=8) :: Z
    end type

    ! Definition of priority queue
    type pqueue
      type(node), allocatable :: buf(:)
      integer :: n = 0
    contains
      procedure :: PQpop
      procedure :: PQpush
      procedure :: shiftdown
    end type

    type (pqueue) :: priorityqueue

contains

    subroutine shiftdown(this, a)
    !*****************************************************************************
    ! This function sort the queue based on increasing elevations priority

      class (pqueue)  :: this
      integer :: a, parent, child

      associate (x => this%buf)
      parent = a

      do while(parent*2 <= this%n)
        child = parent*2
        if (child + 1 <= this%n) then
          if (x(child+1)%Z < x(child)%Z ) then
            child = child +1
          end if
        end if

        if (x(parent)%Z > x(child)%Z) then
          x([child, parent]) = x([parent, child])
          parent = child
        else
          exit
        end if
      end do
      end associate

    end subroutine shiftdown

    function PQpop(this) result (res)
    !*****************************************************************************
    ! This function pops first values in a priority queue

      class(pqueue) :: this
      type(node)   :: res

      res = this%buf(1)
      this%buf(1) = this%buf(this%n)
      this%n = this%n - 1

      call this%shiftdown(1)

    end function PQpop

    subroutine PQpush(this, Z, id)
    !*****************************************************************************
    ! This function pushes new values in a priority queue

      class(pqueue), intent(inout) :: this
      real(kind=8) :: Z
      integer  :: id

      type(node)  :: x
      type(node), allocatable  :: tmp(:)
      integer :: ii

      x%Z = Z
      x%id = id
      this%n = this%n +1

      if (.not.allocated(this%buf)) allocate(this%buf(1))
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if

      this%buf(this%n) = x
      ii = this%n
      do
        ii = ii / 2
        if (ii==0) exit
        call this%shiftdown(ii)
      end do

    end subroutine PQpush

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

    subroutine fillBarnes(elevation, sealevel, demH, pydnodes)

      logical :: change
      integer :: pydnodes, n, k, p
      real(kind=8),intent(in) :: sealevel
      real(kind=8),intent(in) :: elevation(pydnodes)
      real(kind=8),intent(inout) :: demH(pydnodes)


      type (node)  :: pt
      logical :: flag(pydnodes)
      real(kind=8) :: fill(pydnodes)


      ! Push edges to priority queue
      flag = .False.
      fill = elevation
      do k = 1, bds
        call priorityqueue%PQpush(fill(k),k)
        flag(k) = .True.
      enddo

      ! Use priority queue to remove pits
      do while(priorityqueue%n >0)
        pt = priorityqueue%PQpop()
        k = pt%id
        demH(k) = fill(k)
        loop: do p = 1, 20
          if( neighbours(k,p) < 0 ) exit loop
          n = neighbours(k,p)+1
          if(.not. flag(n))then
            flag(n) = .True.
            fill(n) = max(fill(n),fill(k)+eps)
            call priorityqueue%PQpush(fill(n), n)
          endif
        enddo loop
      enddo

      ! Limit elevation based on maximum filling thickness
      do k = 1, pydnodes
        if(elevation(k)>=sealevel)then
          if(demH(k)-elevation(k)>fill_TH) demH(k) = elevation(k) + fill_TH
        else
          ! demH(k) = elevation(k)
          if(demH(k)>sealevel)then
              demH(k) = min(demH(k),sealevel+fill_TH)
          else
            demH(k) = elevation(k)
          endif
        endif
      enddo

      return

    end subroutine fillBarnes

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
