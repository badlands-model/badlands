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
!       Filename:  StratalClass.f90
!
!    Description:  Implements the geometry and topology functions of the stratigraphic grid
!
!        Version:  1.0
!        Created:  02/06/15 07:56:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module stratal_class

  use bilinear
  use topology
  use parallel
  use hydroUtil
  use parameters
  use kdtree2_module
  use kdtree2_precision_module

  implicit none

  ! Stratal grid resolution
  integer::slay,ilay,layNb
  real(kind=8)::time_layer,time_strata

  real(kind=8),dimension(:,:),allocatable::stratadata
  type(kdtree2),pointer::stratatree
  type(kdtree2_result),dimension(10)::strataRslt

  type sediment_parameters
     ! Diameter in millimetres
     real(kind=8)::dia
     ! Density in kg.m-3
     real(kind=8)::dens
     ! Fall velocity
     real(kind=8)::vfall
     ! Marine diffusivity
     real(kind=8)::diffm
     ! Aerial diffusivity
     real(kind=8)::diffa
     ! Erodibility coefficient
     real(kind=8)::ero
  end type sediment_parameters
  type(sediment_parameters),dimension(:),allocatable::sediments

  ! Stratigraphic layer
  character(len=128),dimension(:),allocatable::flayh

  integer,dimension(:),allocatable::lnID
  integer,dimension(:),allocatable::delPt
  integer,dimension(:),allocatable::regPtNb
  integer,dimension(:,:),allocatable::regPtID

  real(kind=8),dimension(:),allocatable::ilayh
  real(kind=8),dimension(:),allocatable::lay_base
  real(kind=8),dimension(:,:),allocatable::lay_thick
  real(kind=8),dimension(:,:),allocatable::lay_porosity
  real(kind=8),dimension(:,:,:),allocatable::lay_sed

  ! Active layer parameters for stratigraphic mesh
  real(kind=8)::active_thick
  real(kind=8),dimension(:),allocatable::alay_thick
  real(kind=8),dimension(:),allocatable::alay_porosity
  real(kind=8),dimension(:,:),allocatable::alay_sed

  ! Active layer parameters for delaunay grid
  real(kind=8),dimension(:,:),allocatable::alay_dsed

  real,dimension(:),allocatable::biX,biY,interpVal
  real,dimension(:,:),allocatable::biZ
  real,dimension(:,:,:,:),allocatable::biSed

  ! Sediment influx
  real(kind=8),dimension(:,:),allocatable::Qs_inS
  real(kind=8),dimension(:,:),allocatable::change_localS

contains

  ! =====================================================================================
  subroutine constructStrata

    logical::found

    integer::i,j,p,n,k,iu

    real(kind=8)::s
    real(kind=8),dimension(totgrn)::sth
    real(kind=8),dimension(upartN,2)::lay_coord
    real(kind=8),dimension(nbnodes,ilay,totgrn)::lay_sed1

    ! Stratal layer dimension
    slay=int((time_end-time_start)/time_layer)+ilay

    ! Read initial deposit
    do k=1,ilay
      if(ilayh(k)==0)then
        inquire(file=flayh(k),exist=found)
        if(.not.found)then
          if(pet_id==0) print*,'Cannot find XmL stratal file ',k
          call mpi_finalize(rc)
        endif
        iu=79
        open(iu,file=flayh(k),status="old",action="read",iostat=rc)
        rewind(iu)
        sth=0.
        do n=1,nbnodes
          read(iu,*)sth(1:totgrn)
          s=sum(sth)
          do p=1,totgrn
            lay_sed1(n,k,p)=sth(p)/s
          enddo
        enddo
      else
        s=ilayh(k)
        do n=1,nbnodes
          lay_sed1(n,k,1:totgrn)=s/totgrn
        enddo
      endif
    enddo

    ! Define coarse regular mesh for interpolation
    if(allocated(biX)) deallocate(biX)
    if(allocated(biY)) deallocate(biY)
    if(allocated(biZ)) deallocate(biZ)
    if(allocated(biSed)) deallocate(biSed)
    allocate(biX(nx+2))
    allocate(biY(ny+2))
    allocate(biZ(nx+2,ny+2))
    allocate(biSed(ilay,totgrn,nx+2,ny+2))

    biX(1)=real(minx-dx)
    do k=2,nx+2
      biX(k)=biX(k-1)+real(dx)
    enddo
    biY(1)=real(miny-dx)
    do k=2,ny+2
      biY(k)=biY(k-1)+real(dx)
    enddo
    p=0
    do j=2,ny+1
      do i=2,nx+1
        p=p+1
        biZ(i,j)=real(coordZ(p))
        do k=1,ilay
          do n=1,totgrn
            biSed(k,n,i,j)=real(lay_sed1(p,k,n))
          enddo
        enddo
      enddo
    enddo
    biZ(1,1)=biZ(2,2)
    biZ(nx+2,1)=biZ(nx+1,2)
    biZ(1,ny+2)=biZ(2,ny+1)

    biZ(nx+2,ny+2)=biZ(nx+1,ny+1)
    biZ(2:nx+1,1)=biZ(2:nx+1,2)
    biZ(2:nx+1,ny+2)=biZ(2:nx+1,ny+1)
    biZ(1,2:nx+1)=biZ(2,2:nx+1)
    biZ(nx+2,2:nx+1)=biZ(nx+1,2:nx+1)

    do k=1,ilay
      do n=1,totgrn
        biSed(k,n,1,1)=biSed(k,n,2,2)
        biSed(k,n,nx+2,1)=biSed(k,n,nx+1,2)
        biSed(k,n,1,ny+2)=biSed(k,n,2,ny+1)
        biSed(k,n,nx+2,ny+2)=biSed(k,n,nx+1,ny+1)
        biSed(k,n,2:nx+1,1)=biSed(k,n,2:nx+1,2)
        biSed(k,n,2:nx+1,ny+2)=biSed(k,n,2:nx+1,ny+1)
        biSed(k,n,1,2:nx+1)=biSed(k,n,2,2:nx+1)
        biSed(k,n,nx+2,2:nx+1)=biSed(k,n,nx+1,2:nx+1)
      enddo
    enddo

    ! Allocate the array dimension
    if(allocated(lay_base)) deallocate(lay_base)
    allocate(lay_base(upartN))
    if(allocated(lay_thick)) deallocate(lay_thick)
    if(allocated(lay_sed)) deallocate(lay_sed)
    allocate(lay_thick(upartN,slay))
    allocate(lay_sed(upartN,slay,totgrn))
    lay_thick=0.0
    lay_sed=0.0

    if(allocated(interpVal)) deallocate(interpVal)
    allocate(interpVal(upartN))

    ! Build grid X-Y coordinates for fine stratal mesh
    do k=1,upartN
      n=unodeID(k)
      lay_coord(k,1)=tcoordX(n)
      lay_coord(k,2)=tcoordY(n)
    enddo

    ! Interpolate values from coarse regular grid to unstructured stratal grid
    call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biZ,upartN,real(lay_coord(1:upartN,1)),real(lay_coord(1:upartN,2)),interpVal)
    lay_base=interpVal

    do k=1,ilay
      do n=1,totgrn
        call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biSed(k,n,1:nx+2,1:ny+2),upartN,real(lay_coord(1:upartN,1)),real(lay_coord(1:upartN,2)),interpVal)
        lay_sed(1:upartN,k,n)=interpVal(1:upartN)
      enddo
    enddo

    do n=1,upartN
      lay_thick(n,1:slay)=0.0
      do k=1,ilay
        s=0.0_8
        do p=1,totgrn
          s=s+lay_sed(n,k,p)
        enddo
        lay_thick(n,k)=s
        lay_base(n)=lay_base(n)-s
      enddo
    enddo

    layNb=ilay

    return

  end subroutine constructStrata
  ! =====================================================================================
end module stratal_class
