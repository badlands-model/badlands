!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements the stack used for the Planchon & Darboux depression filling
! algorithm

module pdstack

  public

  ! Define the data-structure to hold the data
  type stack_var
     integer, allocatable :: data(:)
     integer              :: size = 0
  end type stack_var
  type(stack_var) :: stack1,stack2
      
  ! Set the size of allocated memory blocks
  integer, parameter, private :: block_size = 10000000

contains

  ! Push ----------------------------------------------------------------------
  subroutine push(s, e)
    type(stack_var), intent(inout) :: s
    integer, intent(in)            :: e
    integer, allocatable :: wk(:)
    if (.not. allocated(s%data)) then
       ! Allocate space if not yet done
       allocate(s%data(block_size))

    elseif (s%size == size(s%data)) then
       ! Grow the allocated space
       allocate(wk(size(s%data)+block_size))
       wk(1:s%size) = s%data
       call move_alloc(wk,s%data)

    end if

    ! Store the data in the stack
    s%size = s%size + 1
    s%data(s%size) = e
  end subroutine push

  ! Pop -----------------------------------------------------------------------
  function pop(s)
    integer :: pop
    type(stack_var), intent(inout) :: s
    if (s%size == 0 .or. .not. allocated(s%data)) then
       pop = 0
       return
    end if
    pop = s%data(s%size)
    s%size = s%size - 1
  end function pop

  ! Peek ----------------------------------------------------------------------
  integer function peek(s)
    type(stack_var), intent(inout) :: s
    if (s%size == 0 .or. .not. allocated(s%data)) then
       peek = 0
       return
    end if
    peek = s%data(s%size)
  end function peek

  ! Empty ---------------------------------------------------------------------
  logical function empty(s)
    type(stack_var), intent(inout) :: s
    empty = (s%size == 0 .or. .not. allocated(s%data))
  end function empty

end module pdstack