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
  real(kind=8), allocatable, dimension(:) :: sedb
  real(kind=8), allocatable, dimension(:) :: sedc
  real(kind=8), allocatable, dimension(:) :: sedd
  real(kind=8), allocatable, dimension(:) :: bedb
  real(kind=8), allocatable, dimension(:) :: bedc
  real(kind=8), allocatable, dimension(:) :: bedd

contains


  subroutine erodibility_factor(sedsupply, upslope, sedsupFactor, bedpropor)

      real(kind=8) :: sedsupply, upslope
      real(kind=8) :: sedsupFactor, bedpropor

      if(erofct)then
        ! Interpolate bedload percentage based on upstream slope
        sedsupFactor = ispline(sedsupply, sedsup, sedsupval, sedNb)!, sedb, sedc, sedd, sedNb)
        print*,'sedsup',sedsupply,sedsupFactor
        ! Interpolate erodibility factor from sediment supply
        bedpropor = ispline(upslope, bedslp, bedprop, slpNb)!, bedb, bedc, bedd, slpNb)
        print*,'upslope',upslope,bedpropor
      else
        sedsupFactor = 1.
        bedpropor = 1.
      endif

      return

  end subroutine erodibility_factor

  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  subroutine spline (x, y, b, c, d, n)

      integer :: n
      real(kind=8) :: x(n), y(n), b(n), c(n), d(n)
      integer :: i, j, gap
      real(kind=8) :: h

      gap = n-1
      ! Check input
      if ( n < 2 ) return
      if ( n < 3 ) then
        b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
        c(1) = 0.
        d(1) = 0.
        b(2) = b(1)
        c(2) = 0.
        d(2) = 0.
        return
      endif

      ! Step 1: preparation
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
        d(i) = x(i+1) - x(i)
        b(i) = 2.0*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
      end do

      ! Step 2: end conditions
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if(n /= 3) then
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if

      ! Step 3: forward elimination
      do i = 2, n
        h = d(i-1)/b(i-1)
        b(i) = b(i) - h*d(i-1)
        c(i) = c(i) - h*c(i-1)
      end do

      ! Step 4: back substitution
      c(n) = c(n)/b(n)
      do j = 1, gap
        i = n-j
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do

      ! Step 5: compute spline coefficients
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
      do i = 1, gap
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3.*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)

  end subroutine spline

  ! Function ispline evaluates the cubic spline interpolation at point z
  function ispline(u, x, y, n) !b, c, d, n)

      real(kind=8) :: ispline
      integer :: n
      real(kind=8) :: u, x(n), y(n) !, b(n), c(n), d(n)
      integer :: i, j, k
      real(kind=8) :: dx !, a, bb

      ! If u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
        ispline = y(1)
        return
      end if
      if(u >= x(n)) then
        ispline = y(n)
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
      !print*,i
      !a = (y(i+1)-y(i))/(x(i+1)-x(i))
      !bb = y(i) !- (y(i+1)-y(i)) * x(i) /(x(i+1)-x(i))
      dx = u - x(i)
      ispline = (y(i+1)-y(i))/(x(i+1)-x(i)) * dx + y(i) !y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      !print*,'u',u,a * dx + bb
  end function ispline

end module flowclass
