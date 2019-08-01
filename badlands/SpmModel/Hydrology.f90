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
!       Filename:  Hydrology.f90
!
!    Description:  Compute morphometrics based on hydrological features
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module hydrology

  use sorting
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none

  ! Local arrays
  integer,dimension(:),allocatable::intArray,allocs,donorCount,partStack

  real(kind=8),dimension(:),allocatable::watercell
  real(kind=8),dimension(:),allocatable::nZ,nH,change_local,filldem

contains

  ! =====================================================================================
  subroutine define_landscape_network

    integer,dimension(npets)::disps
    integer,dimension(npets+1)::stpStack
    integer::k,p,j,lowestID,success,partnb,extra

    if(.not.allocated(partStack)) allocate(partStack(npets))
    if(.not.allocated(receivers)) allocate(receivers(dnodes))
    if(.not.allocated(stackOrder)) allocate(stackOrder(dnodes))
    if(.not.allocated(lstackOrder)) allocate(lstackOrder(dnodes))
    if(.not.allocated(donorCount)) allocate(donorCount(dnodes))
    if(.not.allocated(discharge)) allocate(discharge(dnodes))
    if(.not.allocated(bsID)) allocate(bsID(dnodes))
    if(.not.allocated(chi)) allocate(chi(dnodes))

    receivers=-1
    discharge=0.
    watercell=0
    bsID=-1
    chi=0.
    do j=1,upartN
      k=unodeID(j)
      watercell(k)=filldem(k)-spmZ(k)
      receivers(k)=k
      lowestID=k
      do p=1,delaunayVertex(k)%ngbNb
        if(delaunayVertex(k)%ngbID(p)>0)then
          if(filldem(delaunayVertex(k)%ngbID(p))<filldem(lowestID))lowestID=delaunayVertex(k)%ngbID(p)
        endif
      enddo
      receivers(k)=lowestID
      discharge(k)=precipitation(k)*voronoiCell(k)%area
    enddo
    call mpi_allreduce(mpi_in_place,receivers,dnodes,mpi_integer,mpi_max,badlands_world,rc)
    call mpi_allreduce(mpi_in_place,watercell,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
    call mpi_allreduce(mpi_in_place,discharge,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

    baseNb=0
    if(.not.allocated(baselist)) allocate(baselist(dnodes))
    baselist=-1
    donorCount=0
    do k=1,dnodes
       ! Baselevel
       if(receivers(k)==k)then
          baseNb=baseNb+1
          baselist(baseNb)=k
       endif
      donorCount(receivers(k))=donorCount(receivers(k))+1
    enddo

    ! Index of donors number
    if(.not.allocated(indexArray)) allocate(indexArray(dnodes+1))
    indexArray=0
    maxrcvs=0
    indexArray(dnodes+1)=dnodes+1
    do k=dnodes,1,-1
       indexArray(k)=indexArray(k+1)-donorCount(k)
       maxrcvs=max(maxrcvs,indexArray(k+1)-indexArray(k))
    enddo

    ! List of donors
    if(.not.allocated(allocs)) allocate(allocs(dnodes))
    if(.not.allocated(intArray)) allocate(intArray(dnodes))
    if(.not.allocated(donorsList)) allocate(donorsList(dnodes))
    allocs=-1
    intArray=0
    lstackOrder=-1
    do k=1,dnodes
       donorsList(indexArray(receivers(k))+intArray(receivers(k)))=k
       intArray(receivers(k))=intArray(receivers(k))+1
    enddo

    ! Base level partitioning
    stpStack=0
    partStack=0
    if(baseNb<=npets)then
      p=pet_id+1
      partStack(p)=0
      if(p<=baseNb)partStack(p)=1
      stpStack(p+1)=stpStack(p)+partStack(p)
    else
      if(pet_id==0)then
        partnb=int(baseNb/npets)
        extra=mod(baseNb,npets)
        do j=0,npets-1
          partStack(j+1)=partnb
          if(j<extra) partStack(j+1)=partStack(j+1)+1
          stpStack(j+2)=stpStack(j+1)+partStack(j+1)
        enddo
      endif
      call mpi_bcast(stpStack,npets+1,mpi_integer,0,badlands_world,rc)
    endif

    ! Build the ordering stack
    j=0
    do p=stpStack(pet_id+1)+1,stpStack(pet_id+2)
      allocs=-1
      k=baselist(p)
      j=j+1
      lstackOrder(j)=k
      allocs(k)=0
      success=addtostack(p,k,j)
    enddo

    ! Get ordered stack from catchment ID partitioning
    partStack=0
    call mpi_allgather(j,1,mpi_integer,partStack,1,mpi_integer,badlands_world,rc)
    disps=0
    do p=1,npets-1
      disps(p+1)=disps(p)+partStack(p)
    enddo
    call mpi_allgatherv(lstackOrder,j,mpi_integer,stackOrder,partStack,disps,mpi_integer,badlands_world,rc)

    return

  end subroutine define_landscape_network
  ! =====================================================================================
  recursive function addtostack(base,donor,stackID) result(success)

    integer::base,donor,stackID,n,success

    success=1
    do n=indexArray(donor),indexArray(donor+1)-1
      if((allocs(donorsList(n))==-1))then
        donor=donorsList(n)
        stackID=stackID+1
        lstackOrder(stackID)=donor
        allocs(donor)=0
        success=addtostack(base,donor,stackID)
      endif
    enddo

    success=0

  end function addtostack
  ! =====================================================================================
  subroutine planchon_dem_fill_algorithm

    logical::flag
    integer::p,k,l

    real(kind=8)::step

    ! Find the sinks and fill them using Planchon's method
    if(pet_id==0)then
      flag=.true.
      step=1.e-4_8
      do while(flag)
        flag=.false.
        do k=1,dnodes
          if(voronoiCell(k)%border/=1.and.filldem(k)>spmZ(k))then
            ! Loop through the neighbors
            do p=1,delaunayVertex(k)%ngbNb
              l=delaunayVertex(k)%ngbID(p)
              if(l>0)then
                ! In case delaunay elevation greater than neighbor's ones
                if(spmZ(k)>=filldem(l)+step)then
                  filldem(k)=spmZ(k)
                ! Otherwise this is a sink and we perform sink filling
                else
                  if(filldem(k)>filldem(l)+step)then
                    filldem(k)=filldem(l)+step
                    if(filldem(k)-spmZ(k)>fh)then
                      filldem(k)=spmZ(k)+fh
                    else
                      flag=.true.
                    endif
                  endif
                endif
              endif
            enddo
          endif
        enddo
      enddo
    endif
    call mpi_bcast(filldem,dnodes,mpi_double_precision,0,badlands_world,rc)

    return

  end subroutine planchon_dem_fill_algorithm
  ! =====================================================================================
  subroutine compute_vertical_displacement

    integer::k,lid,id

    do lid=1,localNodes
      id=localNodesGID(lid)
      k=stackOrder(id)
      nZ(k)=nZ(k)+tvertDisp(k)*time_step
    enddo

  end subroutine compute_vertical_displacement
  ! =====================================================================================
  subroutine update_grid_borders

    integer::k,p

    if(simulation_time==time_start.or.update3d)then
      if(.not.allocated(delemID)) allocate(delemID(delem))
      delemoo=0
      delemID=-1
      p=0
      do k=1,delem
         if(elemtmask(k)==0.and.uownEID(k)==pet_id)then
           delemoo=delemoo+1
           p=p+1
           delemID(p)=k
         endif
      enddo
    endif

    ! Define the borders according to boundary definition
    filldem=1.e6
    do k=1,dnodes
      ! On border fix DEM elevation
      if(voronoiCell(k)%border==1)then
        spmH(k)=0.0_8
        ! In case there is an outlet
        if(outlet/=0)then
          if(voronoiCell(k)%btype==outlet)then
            spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
          else
             spmZ(k)=1.e6_8
          endif
        else
          !--------- WEST
          ! South West corner
          if(voronoiCell(k)%btype==1)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North West corner
          if(voronoiCell(k)%btype==3)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! West border
          if(voronoiCell(k)%btype==5)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- EAST
          ! South East corner
          if(voronoiCell(k)%btype==2)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North East corner
          if(voronoiCell(k)%btype==4)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! East border
          if(voronoiCell(k)%btype==6)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- SOUTH
          ! South East corner
          if(voronoiCell(k)%btype==2)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! South West corner
          if(voronoiCell(k)%btype==1)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! South border
          if(voronoiCell(k)%btype==7)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- NORTH
          ! North East corner
          if(voronoiCell(k)%btype==4)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North West corner
          if(voronoiCell(k)%btype==3)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North border
          if(voronoiCell(k)%btype==8)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
        endif
        ! Update DEM borders elevation values
        filldem(k)=spmZ(k)
      endif
    enddo

    return

  end subroutine update_grid_borders
  ! =====================================================================================
  subroutine DeriveTrianglePlanes(x,y,id1,id2,id3,z)

    integer::id1,id2,id3
    real(kind=8)::s,x,y,z

    real(kind=8),dimension(3)::d1,d2,n
    real(kind=8),dimension(4)::plane

    d1(1)=tcoordX(id2)-tcoordX(id1)
    d1(2)=tcoordY(id2)-tcoordY(id1)
    d1(3)=tcoordZ(id2)-tcoordZ(id1)

    d2(1)=tcoordX(id3)-tcoordX(id1)
    d2(2)=tcoordY(id3)-tcoordY(id1)
    d2(3)=tcoordZ(id3)-tcoordZ(id1)

    n(1)=d2(2)*d1(3)-d2(3)*d1(2)
    n(2)=d2(3)*d1(1)-d2(1)*d1(3)
    n(3)=d2(1)*d1(2)-d2(2)*d1(1)

    s=1.0/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    plane(1)=n(1)*s
    plane(2)=n(2)*s
    plane(3)=n(3)*s
    plane(4)= -(plane(1)*tcoordX(id1)+plane(2)*tcoordY(id1)+plane(3)*tcoordZ(id1))

    z=-(plane(1)*x+plane(2)*y+plane(4))/plane(3)

  end subroutine DeriveTrianglePlanes
  ! =====================================================================================
  subroutine DeriveTrianglePlanesSed(x,y,id1,id2,id3,sh)

    integer::p,id1,id2,id3
    real(kind=8)::s,x,y,sh,sload1,sload2,sload3

    real(kind=8),dimension(3)::d1,d2,n
    real(kind=8),dimension(4)::plane

    sload1=0.
    sload2=0.
    sload3=0.
    do p=1,flex_lay
      sload1=sload1+ulay_th(id1,p)*(1-ulay_phi(id1,p))*mean_sediment_density+&
        ulay_th(id1,p)*ulay_phi(id1,p)*sea_water_density
      sload2=sload2+ulay_th(id2,p)*(1-ulay_phi(id2,p))*mean_sediment_density+&
        ulay_th(id2,p)*ulay_phi(id2,p)*sea_water_density
      sload3=sload3+ulay_th(id3,p)*(1-ulay_phi(id3,p))*mean_sediment_density+&
        ulay_th(id3,p)*ulay_phi(id3,p)*sea_water_density
    enddo

    ! Sediment thickness
    d1(1)=tcoordX(id2)-tcoordX(id1)
    d1(2)=tcoordY(id2)-tcoordY(id1)
    d1(3)=sload2-sload1

    d2(1)=tcoordX(id3)-tcoordX(id1)
    d2(2)=tcoordY(id3)-tcoordY(id1)
    d2(3)=sload3-sload1

    n(1)=d2(2)*d1(3)-d2(3)*d1(2)
    n(2)=d2(3)*d1(1)-d2(1)*d1(3)
    n(3)=d2(1)*d1(2)-d2(2)*d1(1)

    s=1.0/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    plane(1)=n(1)*s
    plane(2)=n(2)*s
    plane(3)=n(3)*s
    plane(4)= -(plane(1)*tcoordX(id1)+plane(2)*tcoordY(id1)+plane(3)*sload1)

    sh=-(plane(1)*x+plane(2)*y+plane(4))/plane(3)

  end subroutine DeriveTrianglePlanesSed
  ! =====================================================================================
  subroutine DeriveTrianglePlanes2(xy,xa,ya,za,z)

    real(kind=8),dimension(2)::xy
    real(kind=8),dimension(3)::xa,ya,za
    real(kind=8)::s,z

    real(kind=8),dimension(3)::d1,d2,n
    real(kind=8),dimension(4)::plane

    d1(1)=xa(2)-xa(1)
    d1(2)=ya(2)-ya(1)
    d1(3)=za(2)-za(1)

    d2(1)=xa(3)-xa(1)
    d2(2)=ya(3)-ya(1)
    d2(3)=za(3)-za(1)

    n(1)=d2(2)*d1(3)-d2(3)*d1(2)
    n(2)=d2(3)*d1(1)-d2(1)*d1(3)
    n(3)=d2(1)*d1(2)-d2(2)*d1(1)

    s=1.0/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    plane(1)=n(1)*s
    plane(2)=n(2)*s
    plane(3)=n(3)*s
    plane(4)= -(plane(1)*xa(1)+plane(2)*ya(1)+plane(3)*za(1))

    z=-(plane(1)*xy(1)+plane(2)*xy(2)+plane(4))/plane(3)

  end subroutine DeriveTrianglePlanes2
  ! =====================================================================================
  subroutine inside_triangle(xa,ya,xb,yb,l)

    integer,intent(out)::l

    real(kind=8),intent(in)::xa,ya,xb(3),yb(3)
    real(kind=8)::det0,det1,det2

    l=-1

    det0=(xb(2)-xb(1))*(ya-yb(1))-(yb(2)-yb(1))*(xa-xb(1))
    det1=(xb(3)-xb(2))*(ya-yb(2))-(yb(3)-yb(2))*(xa-xb(2))
    det2=(xb(1)-xb(3))*(ya-yb(3))-(yb(1)-yb(3))*(xa-xb(3))

    if(det0>=0.and.det1>=0.and.det2>=0)then
      l=1
    elseif(det0<=0.and.det1<=0.and.det2<=0)then
      l=1
    endif

    return

  end subroutine inside_triangle
  ! =====================================================================================
end module hydrology
