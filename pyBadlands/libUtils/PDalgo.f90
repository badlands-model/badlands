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

  ! Define the data-structure to hold the data
  integer, allocatable :: stack1(:)
  integer, allocatable :: stack2(:)
  integer :: size1 = 0
  integer :: size2 = 0
  integer, parameter, private :: block_size = 10000000

contains

  subroutine push1(e)

      integer, intent(in) :: e
      integer, allocatable :: wk(:)

      if (.not. allocated(stack1)) then
         ! Allocate space if not yet done
         allocate(stack1(block_size))

      elseif (size1 == size(stack1)) then
         ! Grow the allocated space
         allocate(wk(size(stack1)+block_size))
         wk(1:size1) = stack1
         call move_alloc(wk,stack1)
      end if
      ! Store the data in the stack
      size1 = size1 + 1
      stack1(size1) = e

      return

  end subroutine push1

  subroutine push2(e)

      integer, intent(in) :: e
      integer, allocatable :: wk(:)

      if (.not. allocated(stack2)) then
         ! Allocate space if not yet done
         allocate(stack2(block_size))

      elseif (size2 == size(stack2)) then
         ! Grow the allocated space
         allocate(wk(size(stack2)+block_size))
         wk(1:size2) = stack2
         call move_alloc(wk,stack2)
      end if
      ! Store the data in the stack
      size2 = size2 + 1
      stack2(size2) = e

      return

  end subroutine push2

  subroutine filling2(elevation,pyNgbs,fillTH,epsilon,pybounds,sealimit,demH,pydnodes)

      integer :: pydnodes
      integer,intent(in) :: pybounds
      real(kind=8),intent(in) :: sealimit
      real(kind=8),intent(in) :: fillTH
      real(kind=8),intent(in) :: epsilon
      real(kind=8),intent(in) :: elevation(pydnodes)
      integer,intent(in) :: pyNgbs(pydnodes,20)
      real(kind=8),intent(out) :: demH(pydnodes)

      logical :: flag
      integer :: p, k, l
      real(kind=8) :: minH

      size1 = 0
      size2 = 0
      demH(1:pybounds) = elevation(1:pybounds)
      do k=pybounds+1,pydnodes
          if( demH(k) > sealimit) demH(k) = 1.e6
          call push1(k)
      enddo

      flag = .true.
      do while( flag )
          flag = .false.
          do l = 1,size1
              k = stack1(l)
              if( demH(k) > elevation(k) )then
                  minH = demH(pyNgbs(k,1)+1)
                  p = 2
                  do while(pyNgbs(k,p) >= 0)
                      minH = min(demH(pyNgbs(k,p)+1),minH)
                      p = p+1
                  enddo
                  if (elevation(k) >= minH + epsilon )then
                       demH(k) = elevation(k)
                  else
                      if( demH(k) > minH + epsilon )then
                          demH(k) = minH + epsilon
                          if( demH(k) - elevation(k) > fillTH )then
                              demH(k) = elevation(k) + fillTH
                          else
                              call push2(k)
                              flag = .true.
                          endif
                      else
                          call push2(k)
                      endif
                  endif
              endif
          enddo
          if(flag)then
              stack1 = 0
              size1 = size2
              stack1(1:size1) = stack2(1:size2)
              size2 = 0
          endif
       enddo

       return

  end subroutine filling2

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
          if( demH(k) > sealimit)then 
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
