!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements the erodibility function class
module flowclass

  implicit none

  logical :: erofct

  integer :: sedNb
  integer :: slpNb
  real(kind=8), allocatable, dimension(:) :: sedsup
  real(kind=8), allocatable, dimension(:) :: sedsupval
  real(kind=8), allocatable, dimension(:) :: bedslp
  real(kind=8), allocatable, dimension(:) :: bedprop

contains


  subroutine erodibility_factor(sedsupply, upslope, sedsupFactor, bedperc)

      real(kind=8),intent(in) :: sedsupply, upslope
      real(kind=8),intent(out) :: sedsupFactor, bedperc

      if(erofct)then
        ! Interpolate bedload percentage based on upstream slope
        sedsupFactor = linear(sedsupply, sedsup, sedsupval, sedNb)
        ! Interpolate erodibility factor from sediment supply
        bedperc = linear(upslope, bedslp, bedprop, slpNb)
      else
        sedsupFactor = 1.
        bedperc = 1.
      endif

      return

  end subroutine erodibility_factor

  ! Function linear evaluates the linear interpolation at point u
  function linear(u, x, y, n)

      real(kind=8) :: linear
      integer :: n
      real(kind=8) :: u, x(n), y(n)
      integer :: i, j, k

      ! If u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
        linear = y(1)
        return
      end if
      if(u >= x(n)) then
        linear = y(n)
        return
      end if

      ! Binary search for for i, such that x(i) <= u <= x(i+1)
      i = 1
      j = n+1
      do while (j > i+1)
        k = (i+j)/2
        if(u < x(k)) then
          j=k
          else
          i=k
         endif
      enddo

      ! Evaluate linear interpolation
      linear = (y(i+1)-y(i))/(x(i+1)-x(i)) * (u - x(i)) + y(i)

  end function linear

end module flowclass
