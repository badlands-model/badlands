! =====================================================================================
! BADLANDS (BAsin anD LANdscape DynamicS)
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  Bilinear.f90
!
!    Description:  Implements a bilinear interpolotion from regular grid to scattered data
!
!        Version:  1.0
!        Created:  08/05/15 10:49:25
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module bilinear

  implicit none

  real(kind=8),parameter::epsilon=0.1
  real(kind=8),parameter::epsilon2=0.1*0.1

contains

  ! =====================================================================================
  function binarysearch(length,array,value)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "d" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif

    integer,intent(in)::length
    real,dimension(length),intent(in)::array
    real,intent(in)::value

    integer::binarysearch

    integer::left,middle,right
    real::d

    d=1e-9

    left=1
    right=length
    do
      if(left>right)then
        exit
      endif

      middle=nint((left+right)/2.0)
      if(abs(array(middle)-value)<=d)then
        binarySearch=middle
        return
      elseif(array(middle)>value)then
        right=middle-1
      else
        left=middle+1
      endif
    enddo
    binarysearch=right

  end function binarysearch
  ! ============================================================================
  real function interpolate_circulation_grid(x_len,x_array,y_len,y_array,f,x,y)

    integer,intent(in)::x_len,y_len
    real,dimension(x_len),intent(in)::x_array
    real,dimension(y_len),intent(in)::y_array
    real,dimension(x_len,y_len),intent(in)::f
    real,intent(in)::x,y
    real::denom,x1,x2,y1,y2
    integer::i,j

    i=binarysearch(x_len,x_array,x)
    j=binarysearch(y_len,y_array,y)
    x1=x_array(i)
    x2=x_array(i+1)
    y1=y_array(j)
    y2=y_array(j+1)
    denom=(x2-x1)*(y2-y1)
    interpolate_circulation_grid=(f(i,j)*(x2-x)*(y2-y)+f(i+1,j)*(x-x1)*(y2-y)+ &
      f(i,j+1)*(x2-x)*(y-y1)+f(i+1,j+1)*(x-x1)*(y-y1))/denom

  end function interpolate_circulation_grid
  ! =====================================================================================
  subroutine interpolate_grid_bilinear(x_len,x_array,y_len,y_array,f,u_len,x,y,ibv)

    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation

    integer,intent(in)::x_len,y_len,u_len
    real,dimension(x_len),intent(in)::x_array
    real,dimension(y_len),intent(in)::y_array
    real,dimension(x_len,y_len),intent(in)::f
    real,dimension(u_len),intent(in)::x,y
    real,dimension(u_len),intent(out)::ibv

    real::denom,x1,x2,y1,y2
    integer::i,j,p,i2,j2

    do p=1,u_len

      i=binarysearch(x_len,x_array,x(p))
      j=binarysearch(y_len,y_array,y(p))
      i2=i+1
      j2=j+1
      if(i2>x_len)i2=i
      if(j2>y_len)j2=j

      x1=x_array(i)
      x2=x_array(i2)
      y1=y_array(j)
      y2=y_array(j2)
      denom=(x2-x1)*(y2-y1)
      if(denom==0)then
        ibv(p)=f(i,j)
      else
        ibv(p)=(f(i,j)*(x2-x(p))*(y2-y(p))+f(i2,j)*(x(p)-x1)*(y2-y(p))+ &
          f(i,j2)*(x2-x(p))*(y(p)-y1)+f(i2,j2)*(x(p)-x1)*(y(p)-y1))/denom
      endif
    enddo

  end subroutine interpolate_grid_bilinear
  ! =====================================================================================
  subroutine interpolate_grid_bilinear3(x_len,x_array,y_len,y_array,f1,f2,f3,u_len,x,y,ib1,ib2,ib3)

    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation

    integer,intent(in)::x_len,y_len,u_len
    real,dimension(x_len),intent(in)::x_array
    real,dimension(y_len),intent(in)::y_array
    real,dimension(x_len,y_len),intent(in)::f1,f2,f3
    real,dimension(u_len),intent(in)::x,y
    real,dimension(u_len),intent(out)::ib1,ib2,ib3

    real::denom,x1,x2,y1,y2
    integer::i,j,p,i2,j2

    do p=1,u_len

      i=binarysearch(x_len,x_array,x(p))
      j=binarysearch(y_len,y_array,y(p))
      i2=i+1
      j2=j+1
      if(i2>x_len)i2=i
      if(j2>y_len)j2=j

      x1=x_array(i)
      x2=x_array(i2)
      y1=y_array(j)
      y2=y_array(j2)
      denom=(x2-x1)*(y2-y1)
      if(denom==0)then
        ib1(p)=f1(i,j)
        ib2(p)=f2(i,j)
        ib3(p)=f3(i,j)
      else
        ib1(p)=(f1(i,j)*(x2-x(p))*(y2-y(p))+f1(i2,j)*(x(p)-x1)*(y2-y(p))+ &
          f1(i,j2)*(x2-x(p))*(y(p)-y1)+f1(i2,j2)*(x(p)-x1)*(y(p)-y1))/denom
        ib2(p)=(f2(i,j)*(x2-x(p))*(y2-y(p))+f2(i2,j)*(x(p)-x1)*(y2-y(p))+ &
          f2(i,j2)*(x2-x(p))*(y(p)-y1)+f2(i2,j2)*(x(p)-x1)*(y(p)-y1))/denom
        ib3(p)=(f3(i,j)*(x2-x(p))*(y2-y(p))+f3(i2,j)*(x(p)-x1)*(y2-y(p))+ &
          f3(i,j2)*(x2-x(p))*(y(p)-y1)+f3(i2,j2)*(x(p)-x1)*(y(p)-y1))/denom
      endif
    enddo

  end subroutine interpolate_grid_bilinear3
  ! =====================================================================================
  subroutine pointInTriangleBoundingBox(xy,xb,yb,l)

    integer::l
    real(kind=8)::xy(2),xb(3),yb(3)
    real(kind=8)::xmin,ymin,xmax,ymax

    xmin=min(xb(1),min(xb(2),xb(3)))-epsilon
    ymin=min(yb(1),min(yb(2),yb(3)))-epsilon

    xmax=max(xb(1),max(xb(2),xb(3)))+epsilon
    ymax=max(yb(1),max(yb(2),yb(3)))+epsilon

    l=1
    if(xy(1)<xmin.or.xy(2)<ymin.or.xy(1)>xmax.or.xy(2)>ymax) l=0
    return

  end subroutine pointInTriangleBoundingBox
  ! =====================================================================================
  subroutine sideTriangle(xy,x1,y1,x2,y2,dist)

    real(kind=8)::xy(2),x1,y1,x2,y2,dist

    dist=(y2-y1)*(xy(1)-x1)+(-x2+x1)*(xy(2)-y1)

    return

  end subroutine sideTriangle
  ! =====================================================================================
  subroutine naivepointInTriangle(xy,xb,yb,l)

    integer::l
    real(kind=8)::xy(2),xb(3),yb(3),d1,d2,d3

    d1=-1000.
    call sideTriangle(xy,xb(1),yb(1),xb(2),yb(2),d1)
    d2=-1000.
    call sideTriangle(xy,xb(2),yb(2),xb(3),yb(3),d2)
    d3=-1000.
    call sideTriangle(xy,xb(3),yb(3),xb(1),yb(1),d3)

    l=0
    if(d1>=0..and.d2>=0..and.d3>=0.) l=1

    return

  end subroutine naivepointInTriangle
  ! =====================================================================================
  subroutine distanceSquarePt2Segment(xy,x1,y1,x2,y2,dist)

    real(kind=8)::xy(2),x1,y1,x2,y2,dist
    real(kind=8)::squarelgth,squarelgth2,dotProduct

    squarelgth=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
    dotProduct=((xy(1)-x1)*(x2-x1)+(xy(2)-y1)*(y2-y1))/squarelgth

    if(dotProduct<0)then
      dist=(xy(1)-x1)*(xy(1)-x1)+(xy(2)-y1)*(xy(2)-y1)
    elseif(dotProduct<=1)then
      squarelgth2=(x1-xy(1))*(x1-xy(1))+(y1-xy(2))*(y1-xy(2))
      dist=squarelgth2-dotProduct*dotProduct*squarelgth
    else
      dist=(xy(1)-x2)*(xy(1)-x2)+(xy(2)-y2)*(xy(2)-y2)
    endif
    return

  end subroutine distanceSquarePt2Segment
  ! =====================================================================================
  subroutine is_point_in_triangle(xy,xb,yb,l)

    integer,intent(inout)::l

    real(kind=8),intent(in)::xy(2),xb(3),yb(3)
    real(kind=8)::xa,ya,det0,det1,det2

    l=0
    xa=xy(1)
    ya=xy(2)
    det0=(xb(2)-xb(1))*(ya-yb(1))-(yb(2)-yb(1))*(xa-xb(1))
    det1=(xb(3)-xb(2))*(ya-yb(2))-(yb(3)-yb(2))*(xa-xb(2))
    det2=(xb(1)-xb(3))*(ya-yb(3))-(yb(1)-yb(3))*(xa-xb(3))

    if(det0>=0.and.det1>=0.and.det2>=0)then
      l=1
    elseif(det0<=0.and.det1<=0.and.det2<=0)then
      l=1
    endif

    return

  end subroutine is_point_in_triangle
  ! =====================================================================================
  subroutine insideTriangle(xy,xb,yb,l)

    integer::l
    real(kind=8)::xy(2),xb(3),yb(3),dist

    l=0
    call pointInTriangleBoundingBox(xy,xb,yb,l)
    if(l==0) return

    l=0
    call naivepointInTriangle(xy,xb,yb,l)
    if(l==1) return

    dist=10000.
    call distanceSquarePt2Segment(xy,xb(1),yb(1),xb(2),yb(2),dist)
    if(dist<=epsilon2)then
      l=1
      return
    endif

    dist=10000.
    call distanceSquarePt2Segment(xy,xb(2),yb(2),xb(3),yb(3),dist)
    if(dist<=epsilon2)then
      l=1
      return
    endif

    dist=10000.
    call distanceSquarePt2Segment(xy,xb(3),yb(3),xb(1),yb(1),dist)
    if(dist<=epsilon2)then
      l=1
      return
    endif

    l=0
    return

  end subroutine insideTriangle
  ! =====================================================================================
end module bilinear
