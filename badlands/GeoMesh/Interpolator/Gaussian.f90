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
!       Filename:  Gaussian.f90
!
!    Description:  Gaussian smoothing operator to avoid large vertical gradients.
!
!        Version:  1.0
!        Created:  30/06/15 18:25:30
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module gaussian_filter

  use parallel
  use parameters
  use topology

  implicit none

  real(kind=8),dimension(:,:),allocatable::ice_diff,FN,FS,FE,FW,ice_mass

contains

  ! =====================================================================================
  subroutine gaussian_kernel(sigma,kernel)

    real(kind=8),intent(in)::sigma
    real(kind=8),intent(out),dimension(:,:),allocatable::kernel

    real(kind=8),dimension(:,:),allocatable::xl,yl
    integer::radius,i,j
    real(kind=8)::s,trunc

    trunc=4.
    radius=int(trunc*sigma+0.5)
    s=sigma**2

    ! Set up meshgrid.
    if(allocated(kernel)) deallocate(kernel)
    allocate(xl(-radius:radius,-radius:radius))
    allocate(yl(-radius:radius,-radius:radius))
    do j=-radius,radius
      do i=-radius,radius
        xl(i,j)=i
        yl(i,j)=j
      enddo
    enddo

    ! Make kernel.
    allocate(kernel(-radius:radius,-radius:radius))
    kernel=2.0*exp(-0.5*(xl**2+yl**2)/s)
    kernel=kernel/sum(kernel)

    deallocate(xl,yl)

  end subroutine gaussian_kernel
  ! =====================================================================================
  subroutine tile_and_reflect(input,output)
  ! Set up 3x3 tiles around the input.

    real(kind=8),intent(in),dimension(:,:)::input
    real(kind=8),intent(out),dimension(:,:),allocatable::output

    integer::rows,cols

    rows=ubound(input,1)
    cols=ubound(input,2)

    ! Rely on automatic deallocation to clean this up.
    if(allocated(output)) deallocate(output)
    allocate(output(3*rows,3*cols))

    ! There are 3x3 tiles, we start at the top left and set the tiles up row by row.
    ! Top left is flipped left-to-right and up-to-down.
    output(:rows,:cols)=input(rows:1:-1,cols:1:-1)
    ! Top centre is flipped up-to-down
    output(:rows,cols+1:2*cols)=input(rows:1:-1,:)
    ! Top right is flipped left-to-right and up-to-down.
    output(:rows,2*cols+1:3*cols)=input(rows:1:-1,cols:1:-1)
    ! Middle left flipped left-to-right
    output(rows+1:2*rows,:cols)=input(:,cols:1:-1)
    ! Middle centre unchanged
    output(rows+1:2*rows,cols+1:2*cols)=input(:,:)
    ! Middle right flipped left-to-right
    output(rows+1:2*rows,2*cols+1:3*cols)=input(:,cols:1:-1)
    ! Bottom left flipped left-to-right and up-to-down
    output(2*rows+1:3*rows,:cols)=input(rows:1:-1,cols:1:-1)
    ! Bottom cente flipped up-to-down
    output(2*rows+1:3*rows,cols+1:2*cols)=input(rows:1:-1,:)
    ! Bottom right flipped left-to-right and up-to-down
    output(2*rows+1:3*rows,2*cols+1:3*cols)=input(rows:1:-1,cols:1:-1)

  end subroutine tile_and_reflect
  ! =====================================================================================
  subroutine convolve(input,weights,output)
  ! Convolution.

    real(kind=8),intent(in),dimension(:,:)::input,weights
    real(kind=8),intent(inout),dimension(:,:)::output

    ! These are allocated within tile_and_reflect, we rely on automatic
    ! deallocation at the end of the subroutine.
    real(kind=8),dimension(:,:),allocatable,target::tiled_input
    real(kind=8),dimension(:,:),pointer::overlapping

    integer::rows,cols,hw_row,hw_col,i,j,tj,ti

    ! First step is to tile the input.
    rows=ubound(input,1)
    cols=ubound(input,2)
    ! Stands for half weights row.
    hw_row=ubound(weights,1)/2
    hw_col=ubound(weights,2)/2

    ! Only one reflection is done on each side so the weights array cannot be
    ! bigger than width/height of input +1.
    call assert(ubound(weights,1)<rows+1,'Input size too small for weights matrix')
    call assert(ubound(weights,2)<cols+1,'Input size too small for weights matrix')

    ! This ensures that in the masked case, all masked points remain unchanged.
    output(:,:)=input(:,:)
    call tile_and_reflect(input,tiled_input)

    do j=1,cols
      do i=1,rows
        ! Use i, j to offset into equivalent part of the tiled arrays.
        ti=i+rows
        tj=j+cols
        overlapping=>tiled_input(ti-hw_row:ti+hw_row,tj-hw_col:tj+hw_col)
        output(i,j)=sum(weights(:,:)*overlapping)
      enddo
    enddo

    if(allocated(tiled_input))deallocate(tiled_input)

  end subroutine convolve
  ! =====================================================================================
  subroutine assert(statement,msg)

    logical,intent(in)::statement
    character(len=*),intent(in)::msg

    if(.not.statement)then
      write(6,*) msg
      stop 'Assert triggered, see stderr.'
    endif

  end subroutine assert
  ! =====================================================================================
  subroutine gaussian_operator(input,initialkernel,output)

    logical::initialkernel

    real(kind=8)::sigma
    real(kind=8), dimension(:,:),allocatable::kernel
    real(kind=8),intent(in),dimension(:,:)::input
    real(kind=8),intent(inout),dimension(:,:)::output

    if(initialkernel)then
      sigma=0.3*(gausskernelSize*0.5-1)+0.8
      call gaussian_kernel(sigma,kernel)
    endif

    call convolve(input,kernel,output)

  end subroutine gaussian_operator
  ! =====================================================================================
  subroutine ICE_grid

    integer::step,i,j,ic,jc,p,m
    real(kind=8),dimension(nbnodes)::tempZ

    if(ice_dx<dx) ice_dx=dx
    step=int(ice_dx/dx)
    nbix=int(nx/step+1)+2
    nbiy=int(ny/step+1)+2

    ! South west corner
    ice_xo=minx-ice_dx
    ice_yo=miny-ice_dx

    if(ice_xo+(nbix-1)*ice_dx>minx+nx*dx)nbix=nbix-1
    if(ice_yo+(nbiy-1)*ice_dx>miny+ny*dx)nbiy=nbiy-1

    if(allocated(iceX)) deallocate(iceX)
    allocate(iceX(nbix))
    do i=1,nbix
      iceX(i)=(i-1)*ice_dx+ice_xo
    enddo

    if(allocated(iceY)) deallocate(iceY)
    allocate(iceY(nbiy))
    do i=1,nbiy
      iceY(i)=(i-1)*ice_dx+ice_yo
    enddo

    if(allocated(iceZ)) deallocate(iceZ)
    allocate(iceZ(nbix,nbiy))

    if(allocated(iceZb)) deallocate(iceZb)
    allocate(iceZb(nbix,nbiy))

    if(allocated(iceH)) deallocate(iceH)
    allocate(iceH(nbix,nbiy))
    iceH=0.0

    if(allocated(iceU)) deallocate(iceU)
    allocate(iceU(nbix,nbiy))
    iceU=0.0

    ! Read changes in topographic regular grid.
    p=0
    m=0
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
          tempZ(m)=rcoordZ(p)
        endif
      enddo
    enddo

    j=1
    do jc=2,nbiy-1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbix-1
        iceZ(ic,jc)=tempZ(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    iceZ(2:nbix-1,1)=iceZ(2:nbix-1,2)+(iceZ(2:nbix-1,2)-iceZ(2:nbix-1,3))
    iceZ(2:nbix-1,nbiy)=iceZ(2:nbix-1,nbiy-1)+(iceZ(2:nbix-1,nbiy-1)-iceZ(2:nbix-1,nbiy-2))
    iceZ(1,2:nbiy-1)=iceZ(2,2:nbiy-1)+iceZ(2,2:nbiy-1)-iceZ(3,2:nbiy-1)
    iceZ(nbix,2:nbiy-1)=iceZ(nbix-1,2:nbiy-1)+iceZ(nbix-1,2:nbiy-1)-iceZ(nbix-2,2:nbiy-1)

    ! Update corner
    iceZ(1,1)=iceZ(2,2)+iceZ(2,2)-iceZ(3,3)
    iceZ(1,nbiy)=iceZ(2,nbiy-1)+iceZ(2,nbiy-1)-iceZ(3,nbiy-2)
    iceZ(nbix,1)=iceZ(nbix-1,2)+iceZ(nbix-1,2)-iceZ(nbix-2,3)
    iceZ(nbix,nbiy)=iceZ(nbix-1,nbiy-1)+iceZ(nbix-1,nbiy-1)-iceZ(nbix-2,nbiy-2)

  end subroutine ICE_grid
  ! =====================================================================================
  subroutine FLEX_grid

    integer::step,i

    if(flex_dx<dx) flex_dx=dx
    step=int(flex_dx/dx)
    nbfx=int(nx/step)
    nbfy=int(ny/step)

    ! South west corner
    flex_xo=minx-flex_dx
    flex_yo=miny-flex_dx

    if(flex_xo+(nbfx+2)*flex_dx<minx+nx*dx)nbfx=nbfx+1
    if(flex_yo+(nbfy+2)*flex_dx<miny+ny*dx)nbfy=nbfy+1

    if(allocated(flexX)) deallocate(flexX)
    allocate(flexX(nbfx+2))
    do i=1,nbfx+2
      flexX(i)=(i-1)*flex_dx+flex_xo
    enddo

    if(allocated(flexY)) deallocate(flexY)
    allocate(flexY(nbfy+2))
    do i=1,nbfy+2
      flexY(i)=(i-1)*flex_dx+flex_yo
    enddo

    if(allocated(flexZ)) deallocate(flexZ)
    allocate(flexZ(nbfx+2,nbfy+2))

    if(allocated(flexSed)) deallocate(flexSed)
    allocate(flexSed(nbfx+2,nbfy+2))

    if(allocated(load)) deallocate(load)
    allocate(load(nbfx+2,nbfy+2))
    load=0.

    if(allocated(prevload)) deallocate(prevload)
    allocate(prevload(nbfx+2,nbfy+2))

    if(allocated(flexDisp)) deallocate(flexDisp)
    allocate(flexDisp(nbfx+2,nbfy+2))

    ! Declare parameters
    sea_water_density=1025.
    cst1=9.81/flex_rigid
    cst2=cst1*(mean_mantle_density-sea_water_density)
    cst3=cst1*mean_mantle_density
    torb=10.0**(-3.0+2.0_8*(23.0_8-log10(flex_rigid))/3.0_8)  &
      *(dx/27000.0_8)**2

  end subroutine FLEX_grid
  ! =====================================================================================
end module gaussian_filter
! =====================================================================================
