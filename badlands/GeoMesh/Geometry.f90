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
!       Filename:  Geometry.f90
!
!    Description:  Implements the geometry and topology functions
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module geomesh

  use parallel
  use parameters
  use readtopo
  use external_forces
  use readforces
  use readgeomorpho
  use topology
  use outsurf
  use meshPartition
  use m_mrgrnk
  use orderpack
  use restart

  implicit none

contains
  ! =====================================================================================
  subroutine GeoMesher

    logical::found

!     real(kind=8)::time1,time2
    real_type=mpi_real
    max_type=mpi_max
    totgrn=0
    update3d=.false.
    inquire(file=xmlfile,exist=found)
    if(.not.found)then
       if(pet_id==0) print*,'Cannot find XmL input file'
       call mpi_finalize(rc)
    endif

    ! GeoMesh main utilities
!     time1=mpi_wtime()
    call topology_parser
    call geomorphology_parser
    call forces_parser
    call ReadRegular
    if(pet_id==0)then
      if(disp3d.and.restartFlag)then
        call DelaunayTransformForce
      else
        call DelaunayTransform
      endif
    endif
    call mpi_barrier(badlands_world,rc)
    call ReadTriangle
    call DelaunayVoronoiDuality
    call DelaunayBorders
!     if(pet_id==0)then
!       call delaunay_xmf
!       call voronoi_xmf
!     endif

    ! Unstructured Grid Partitioning
    call UnstructureGridPart
    ! Structured Grid Partitioning
    call StructureGridPart

!     time2=mpi_wtime()
!     if(pet_id==0) print*,'Number of points',dnodes
!     if(pet_id==0) print*,'BADLANDS Meshes Declaration & Partitioning Completed (s) ',time2-time1
    simulation_time=time_start

    return

  end subroutine GeoMesher
  ! =====================================================================================
  subroutine ReadRegular

    logical::first

    character(len=1)::line

    integer::iu,id,k,p,m,nnx,nny
    integer::n
    integer,parameter::maxrecs=1000000000

    real(kind=8)::xo,yo

    iu=97
    open(iu,file=regularfile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open regular input file'
      call mpi_finalize(rc)
    endif

    ! Determine total number of lines in file
    nbnodes=0
    nodes_count: do n=1,maxrecs
       read(iu,*,iostat=rc) line
       if(rc/=0) exit nodes_count
       if(n==maxrecs)then
          print*,'Maximum number of nodes number exceeded in the regular input file'
          call mpi_finalize(rc)
       endif
       nbnodes=nbnodes+1
    enddo nodes_count
    rewind(iu)

    ! Allocate node coords
    if(allocated(coordX))deallocate(coordX)
    if(allocated(coordY))deallocate(coordY)
    if(allocated(coordZ))deallocate(coordZ)
    allocate(coordX(nbnodes),coordY(nbnodes),coordZ(nbnodes))
    nx=0
    ny=0
    first=.true.
    do n=1,nbnodes
       read(iu,*)coordX(n),coordY(n),coordZ(n)
       if(n==1)then
          minx=coordX(n)
          miny=coordY(n)
       elseif(n==nbnodes)then
          maxx=coordX(n)
          maxy=coordY(n)
       endif
       if(n==2)dx=coordX(n)-minx
       if(first.and.coordY(n)>miny)then
          nx=n-1
          ny=nbnodes/nx
          first=.false.
       endif
    enddo

    close(iu)

    if(del_area<dx**2.)then
      if(pet_id==0)then
        print*,'ERROR: the chosen XmL delaunay area element is too small.'
        print*,'The delaunay area has been set to',del_area
        print*,'The imported grid has a cell area of',dx**2.0
        print*,'The delaunay area should not be greater than the imported cell area.'
        print*,'If you wish to increase the simulation area in some places define some refine area'
      endif
      call mpi_finalize(rc)
    endif

    ! Create a halo around the number of regular grid
    nnx=nx+2
    nny=ny+2
    bnbnodes=nnx*nny

    ! Define regular grid with additional halo
    if(allocated(rcoordX))deallocate(rcoordX)
    if(allocated(rcoordY))deallocate(rcoordY)
    if(allocated(rcoordZ))deallocate(rcoordZ)
    if(allocated(rtectoZ))deallocate(rtectoZ)
    if(allocated(rainVal))deallocate(rainVal)
    if(allocated(rDisp))deallocate(rDisp)
    if(allocated(rvertDisp))deallocate(rvertDisp)
    if(allocated(rsedload))deallocate(rsedload)
    allocate(rcoordX(bnbnodes),rcoordY(bnbnodes),rcoordZ(bnbnodes),rainVal(bnbnodes))
    allocate(rsedload(bnbnodes))
    if(.not.disp3d)then
      allocate(rtectoZ(bnbnodes),rDisp(bnbnodes,1),rvertDisp(bnbnodes))
    else
      if(allocated(rhyDisp))deallocate(rhyDisp)
      if(allocated(rhxDisp))deallocate(rhxDisp)
      allocate(rhyDisp(bnbnodes),rhxDisp(bnbnodes),rvertDisp(bnbnodes))
      allocate(rtectoZ(bnbnodes),rDisp(bnbnodes,3))
    endif
    xo=minx-dx
    yo=miny-dx
    m=0
    id=0
    rDisp=0.0
    rvertDisp=0.0
    if(rain_event==0)then
      rainVal=1.0
    else
      rainVal=0.0
    endif
    do k=1,nny
      do p=1,nnx
        m=m+1
        rcoordX(m)=xo+(p-1)*dx
        rcoordY(m)=yo+(k-1)*dx
        if(p>1.and.p<nnx.and.k>1.and.k<nny)then
          id=id+1
          rcoordZ(m)=coordZ(id)
        endif
      enddo
    enddo

    ! Get the regular grid corners values
    rcoordZ(1)=coordZ(1)
    rcoordZ(nnx)=coordZ(nx)
    rcoordZ(bnbnodes-nnx+1)=coordZ(nbnodes-nx+1)
    rcoordZ(bnbnodes)=coordZ(nbnodes)

    ! Get the South/North topographic values
    id=0
    do k=1,nny
      do p=1,nnx
        id=id+1
        if(k==1.and.p>1.and.p<nnx)then
          rcoordZ(id)=rcoordZ(id+nnx)
        elseif(k==nny.and.p>1.and.p<nnx)then
          rcoordZ(id)=rcoordZ(id-nnx)
        endif
      enddo
    enddo

    ! Get the East/West topographic values
    id=0
    do k=1,nny
      do p=1,nnx
        id=id+1
        if(p==1.and.k>1.and.k<nny)then
          rcoordZ(id)=rcoordZ(id+1)
        elseif(p==nnx.and.k>1.and.k<nny)then
          rcoordZ(id)=rcoordZ(id-1)
        endif
      enddo
    enddo

    rtectoZ=rcoordZ

    ! Allocate the elements ID for the regular grid
    id = 1
    k = 1
    nbelm=(nnx-1)*(nny-1)
    if(allocated(relmt))deallocate(relmt)
    allocate(relmt(nbelm,4))
    do p=1,nbelm
      relmt(p,1)=id
      relmt(p,2)=id+1
      relmt(p,3)=nnx+id+1
      relmt(p,4)=nnx+id
      id=id+1
      if(id==k*nnx)then
        id=k*nnx+1
        k=k+1
      endif
    enddo

    return

  end subroutine ReadRegular
  ! =====================================================================================
  subroutine DelaunayTransform

    character(len=128)::TINCfile,OPTC,stg

    integer::iu,n,p,k,tnb,issue
    integer,dimension(:,:),allocatable::ptid

    iu=65
    TINfile='TIN.poly'
    call addpath2(TINfile)
    TINCfile='TIN.poly'
    call addpath2(TINCfile)
    TINCfile(len(TINfile):len(TINfile))=CHAR(0)
    open(iu,file=TINfile,status="replace",action="write",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open Triangle Polygon Shape File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    ! A box with four points in 2D, no attributes, one boundary marker.
    tnb=4+2*nx+2*(ny-2)
    if(refineNb>0) tnb=tnb+4*refineNb
    write(iu,'(I10,1X,3(I2,1X))') tnb,2,0,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 1,minx-dx,miny-dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 2,minx-dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 3,maxx+dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 4,maxx+dx,miny-dx,1

    ! Write simulation area boundary
    p=4
    do n=1,nx
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,miny,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,maxy,0
    enddo
    do n=2,ny-1
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx,miny+(n-1)*dx,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,maxx,miny+(n-1)*dx,0
    enddo

    ! Box subject to refinement
    if(refineNb>0)then
      issue=0
      if(allocated(ptid)) deallocate(ptid)
      allocate(ptid(refineNb,4))
      ptid=-1
      k=1
      do n=1,refineNb
        p=p+1
        write(iu,'(I10,1X,2(F16.3,1X),I2)')p,refine_grid(n)%xcoord(1),refine_grid(n)%ycoord(1),0
        ptid(n,k)=p
        k=k+1
        p=p+1
        write(iu,'(I10,1X,2(F16.3,1X),I2)')p,refine_grid(n)%xcoord(1),refine_grid(n)%ycoord(2),0
        ptid(n,k) = p
        k=k+1
        p=p+1
        write(iu,'(I10,1X,2(F16.3,1X),I2)')p,refine_grid(n)%xcoord(2),refine_grid(n)%ycoord(2),0
        ptid(n,k) = p
        k=k+1
        p=p+1
        write(iu,'(I10,1X,2(F16.3,1X),I2)')p,refine_grid(n)%xcoord(2),refine_grid(n)%ycoord(1),0
        ptid(n,k) = p
        k=k+1
        if(refine_grid(n)%xcoord(1)<=minx.or.refine_grid(n)%xcoord(2)>=maxx) issue=1
        if(refine_grid(n)%ycoord(1)<=miny.or.refine_grid(n)%ycoord(2)>=maxy) issue=1
        if(issue==1)then
          if(pet_id==0)then
            print*,'ERROR: one or more the borders/corners of the chosen refined area is not properly set.'
            print*,'The border/corner of the refined area region',k
            print*,'should not be along the border of the regular imported grid.'
            print*,refine_grid(n)%xcoord(1:2),minx,maxx
            print*,refine_grid(n)%ycoord(1:2),miny,maxy
          endif
          call mpi_finalize(rc)
        endif
      enddo
    endif

    ! Segments of the outer box
    tnb=4
    if(refineNb>0) tnb=tnb+4*refineNb
    write(iu,'(I10,1X,I2)') tnb,1
    do n=1,4
       if(n<4)then
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,n+1,1
       else
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,1,1
       endif
    enddo

    ! Segments of the inner refined regions
    if(refineNb>0)then
      do k=1,refineNb
        n=n+1
        write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,ptid(k,1),ptid(k,2),0
        n=n+1
        write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n, ptid(k,2),ptid(k,3),0
        n=n+1
        write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,ptid(k,3),ptid(k,4),0
        n=n+1
        write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,ptid(k,4),ptid(k,1),0
      enddo
    endif

    ! Hole not present
    write(iu,'(I10)') 0

    ! Refinement region area constrain
    if(refineNb>0)then
      write(iu,'(I10)') refineNb
      do k=1,refineNb
        write(iu,'(I10,1X,F16.3,1X,F16.3,1X,F16.3)')k,refine_grid(k)%xcoord(1)+10.0_8,&
          refine_grid(k)%ycoord(1)+10.0_8,refine_grid(k)%area
      enddo
    endif
    close(iu)

    ! Conforming Delaunay Triangulation options
    OPTC="-q"
    call append_nbreal(OPTC,del_angle)
    if(refineNb==0)then
      stg="a"
    else
      stg="aa"
    endif
    call append_str(OPTC,stg)
    call append_nbreal(OPTC,del_area)
    OPTC(len(OPTC):len(OPTC))=CHAR(0)

    ! C-function Triangle call (Shewchuck's algorithm)
    call trianglegen(OPTC,TINCfile)

    return

  end subroutine DelaunayTransform
  ! =====================================================================================
  subroutine DelaunayTransformForce

    character(len=128)::TINCfile,OPTC,stg

    integer::iu,n,p,tnb,l

    call getSPM_hdf5topography3D

    iu=65
    TINfile='TIN.poly'
    call addpath2(TINfile)
    TINCfile='TIN.poly'
    call addpath2(TINCfile)
    TINCfile(len(TINfile):len(TINfile))=CHAR(0)
    open(iu,file=TINfile,status="replace",action="write",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open Triangle Polygon Shape File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    ! A box with four points in 2D, no attributes, one boundary marker.
    tnb=4+2*nx+2*(ny-2)+newnds
    write(iu,'(I10,1X,3(I2,1X))') tnb,2,0,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 1,minx-dx,miny-dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 2,minx-dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 3,maxx+dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 4,maxx+dx,miny-dx,1

    ! Write simulation area boundary
    p=4
    do n=1,nx
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,miny,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,maxy,0
    enddo
    do n=2,ny-1
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx,miny+(n-1)*dx,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,maxx,miny+(n-1)*dx,0
    enddo

    do l=1,newnds
      p=p+1
      write(iu,'(I10,1X,2(F16.3,1X),I2)') p,newcoord(l,1:2),0
    enddo

    ! Segments of the outer box
    tnb=4
    write(iu,'(I10,1X,I2)') tnb,1
    do n=1,4
       if(n<4)then
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,n+1,1
       else
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,1,1
       endif
    enddo

    ! Hole not present
    write(iu,'(I10)') 0
    close(iu)

    ! Conforming Delaunay Triangulation options
    OPTC="-q"
    call append_nbreal(OPTC,del_angle)
    stg="a"
    call append_str(OPTC,stg)
    call append_nbreal(OPTC,del_area)
    OPTC(len(OPTC):len(OPTC))=CHAR(0)

    ! C-function Triangle call (Shewchuck's algorithm)
    call trianglegen(OPTC,TINCfile)

    return

  end subroutine DelaunayTransformForce
  ! =====================================================================================
  subroutine ReadTriangle

    integer::iu
    integer::n,a,b,c,m,id

    character(len=256)::trifile

    real(kind=8),dimension(2)::center
    real(kind=8),dimension(:,:),allocatable::map

    ! Read Delaunay nodes
    iu=45
    trifile='TIN.1.node'
    call addpath2(trifile)
    open(iu,file=trifile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open namelist file Triangle Node File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    read(iu,*) dnodes,a,b,c
    if(allocated(tbound)) deallocate(tbound)
    if(allocated(tcoordX)) deallocate(tcoordX)
    if(allocated(tcoordY)) deallocate(tcoordY)
    if(allocated(tcoordZ)) deallocate(tcoordZ)
    if(allocated(tvertDisp)) deallocate(tvertDisp)
    if(allocated(tflex)) deallocate(tflex)
    if(flexure.and.allocated(ulay_th)) deallocate(ulay_th)
    if(flexure.and.allocated(ulay_phi)) deallocate(ulay_phi)
    if(allocated(gtflex)) deallocate(gtflex)
    if(allocated(sedthick)) deallocate(sedthick)
    if(allocated(ice_H)) deallocate(ice_H)
    if(allocated(ice_V)) deallocate(ice_V)
    if(allocated(precipitation)) deallocate(precipitation)
    allocate(tcoordX(dnodes),tcoordY(dnodes),tcoordZ(dnodes),tbound(dnodes))
    allocate(precipitation(dnodes),tvertDisp(dnodes),ice_V(dnodes),ice_H(dnodes))
    allocate(tflex(dnodes),gtflex(dnodes),sedthick(dnodes))
    if(flexure.and..not.restartFlag) allocate(ulay_th(dnodes,flex_lay),ulay_phi(dnodes,flex_lay))
    tcoordZ=0.0
    tflex=0.0
    gtflex=0.0
    tvertDisp=0.0
    ice_H=0.0
    ice_V=0.0
    if(flexure.and..not.restartFlag) ulay_th=0.0
    if(flexure.and..not.restartFlag) ulay_phi=0.0
    if(rain_event==0) precipitation=1.0
    allocate(delaunayVertex(dnodes))
    do n=1,dnodes
       read(iu,*) a,tcoordX(n),tcoordY(n),tbound(n)
    enddo
    close(iu)

    ! Read Delaunay elements
    trifile='TIN.1.ele'
    call addpath2(trifile)
    open(iu,file=trifile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open namelist file Triangle Element File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    read(iu,*) delem,a,b
    if(allocated(delmt)) deallocate(delmt)
    if(allocated(dcentroid)) deallocate(dcentroid)
    if(allocated(Fdata)) deallocate(Fdata)
    if(allocated(elemtmask)) deallocate(elemtmask)
    allocate(delmt(delem,3),dcentroid(delem,2),elemtmask(delem))
    allocate(Fdata(2,delem))
    delemo=0
    elemtmask=0

    ! Find border delaunay elements for output
    do n=1,delem
       read(iu,*) a,delmt(n,1:3)
       center=cmp_centroid(n)
       dcentroid(n,1)=center(1)
       dcentroid(n,2)=center(2)
       do m=1,3
          id=delmt(n,m)
          if(tcoordX(id)<minx.or.tcoordX(id)>maxx) elemtmask(n)=1
          if(tcoordY(id)<miny.or.tcoordY(id)>maxy) elemtmask(n)=1
       enddo
       Fdata(1,n)=dcentroid(n,1)
       Fdata(2,n)=dcentroid(n,2)
       if(elemtmask(n)==0) delemo=delemo+1
    enddo
    close(iu)

    ! Read Delaunay edges
    trifile='TIN.1.edge'
    call addpath2(trifile)
    open(iu,file=trifile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open namelist file Triangle Edge File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    read(iu,*) dedge,a
    allocate(dedg(dedge,3))

    do n=1,dedge
       read(iu,*) a,dedg(n,1:3)
    enddo
    close(iu)

    ! Read Voronoi nodes
    trifile='TIN.1.v.node'
    call addpath2(trifile)
    open(iu,file=trifile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open namelist file Voronoi Nodes File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    read(iu,*) vnodes,a,b,c
    if(allocated(vcoordX)) deallocate(vcoordX)
    if(allocated(vcoordY)) deallocate(vcoordY)
    if(allocated(mask)) deallocate(mask)
    allocate(vcoordX(vnodes),vcoordY(vnodes),mask(vnodes))

    do n=1,vnodes
       read(iu,*) mask(n),vcoordX(n),vcoordY(n)
    enddo
    close(iu)

    ! Mark any duplicate points
    if(allocated(map)) deallocate(map)
    allocate(map(2,vnodes))
    map(1,:)=anint(vcoordX(:)*100000.0)/100000.0
    map(2,:)=anint(vcoordY(:)*100000.0)/100000.0
    call findmap(map,mask)

    ! Read Voronoi edges
    trifile='TIN.1.v.edge'
    call addpath2(trifile)
    open(iu,file=trifile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open namelist file Voronoi Edges File'
      call mpi_finalize(rc)
    endif
    rewind(iu)

    read(iu,*) vedge,a
    allocate(vedg(vedge,2))
    do n=1,vedge
       read(iu,*)a,vedg(n,1:2)
       if(vedg(n,1)>0) vedg(n,1)=mask(vedg(n,1))
       if(vedg(n,2)>0) vedg(n,2)=mask(vedg(n,2))
    enddo
    close(iu)

    return

  end subroutine ReadTriangle
  ! =====================================================================================
  subroutine DelaunayVoronoiDuality

    logical::rec,rec2

    integer::cell,n,id,id1,id2,id3,p,k,nvert,m,tvert,pt(2),tt
    integer,dimension(20)::vsort,sortedID,vorPtsn,vorPts,tsort

    integer,dimension(:),allocatable::ed1,ed2,rk1,rk2

    real(kind=8)::dist,area
    real(kind=8),dimension(20)::vnx,vny,tnx,tny

    ! Delaunay and voronoi parametrisation
    allocate(voronoiCell(dnodes))
    vcellIN=0
    velemIN=0

    ! Rank the edges by point IDs
    if(allocated(ed1)) deallocate(ed1,ed2)
    if(allocated(rk1)) deallocate(rk1,rk2)
    allocate(ed1(dedge),ed2(dedge))
    allocate(rk1(dedge),rk2(dedge))
    ed1(:)=dedg(:,1)
    ed2(:)=dedg(:,2)
    call mrgrnk(ed1,rk1)
    call mrgrnk(ed2,rk2)
    id1=0
    id2=0
    do cell=1,dnodes

      delaunayVertex(cell)%ngbNb=0
      delaunayVertex(cell)%ngbID=-1
      delaunayVertex(cell)%voronoi_edge=0.0_8
      voronoiCell(cell)%perimeter=0.0_8
      voronoiCell(cell)%vertexID=-1
      voronoiCell(cell)%vertexNb=0
      voronoiCell(cell)%border=0
      voronoiCell(cell)%area=0.0_8

      lp_ed1:do n=id1+1,dedge
        if(dedg(rk1(n),1)==cell)then
          id1=n
          rec=.true.
          if(dedg(rk1(n),2)>0)then
            ! Check that the points has not been recorded yet
            do id=1,delaunayVertex(cell)%ngbNb
              if(delaunayVertex(cell)%ngbID(id)==dedg(rk1(n),2)) rec=.false.
            enddo
            if(rec)then
              id=delaunayVertex(cell)%ngbNb+1
              delaunayVertex(cell)%ngbNb=id
              delaunayVertex(cell)%ngbID(id)=dedg(rk1(n),2)
            endif
          endif

          ! Get the Voronoi cell by duality
          rec=.true.
          rec2=.true.
          if(voronoiCell(cell)%vertexNb>=1)then
            do p=1,voronoiCell(cell)%vertexNb
              if(voronoiCell(cell)%vertexID(p)==vedg(rk1(n),1))then
                rec=.false.
              endif
              if(voronoiCell(cell)%vertexID(p)==vedg(rk1(n),2))then
                rec2=.false.
              endif
            enddo
            if(rec.and.vedg(rk1(n),1)>0)then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk1(n),1)
            endif
            if(rec2.and.vedg(rk1(n),2)>0.and.vedg(rk1(n),2)/=vedg(rk1(n),1))then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk1(n),2)
            endif
          else
            if(vedg(rk1(n),1)>0)then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk1(n),1)
            endif
            if(vedg(rk1(n),2)>0.and.vedg(rk1(n),2)/=vedg(rk1(n),1))then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk1(n),2)
            endif
          endif
          if(vedg(rk1(n),1)<=0) voronoiCell(cell)%border=1
        else
          exit lp_ed1
        endif
      enddo lp_ed1

      lp_ed2:do n=id2+1,dedge
        if(dedg(rk2(n),2)==cell)then
          id2=n
          rec=.true.
          if(dedg(rk2(n),1)>0)then
            ! Check that the points has not been recorded yet
            do id=1,delaunayVertex(cell)%ngbNb
              if(delaunayVertex(cell)%ngbID(id)==dedg(rk2(n),1)) rec=.false.
            enddo
            if(rec)then
              id=delaunayVertex(cell)%ngbNb+1
              delaunayVertex(cell)%ngbNb=id
              delaunayVertex(cell)%ngbID(id)=dedg(rk2(n),1)
            endif
          endif

          ! Get the Voronoi cell by duality
          rec=.true.
          rec2=.true.
          if(voronoiCell(cell)%vertexNb>=1)then
            do p=1,voronoiCell(cell)%vertexNb
              if(voronoiCell(cell)%vertexID(p)==vedg(rk2(n),1))then
                rec=.false.
              endif
              if(voronoiCell(cell)%vertexID(p)==vedg(rk2(n),2))then
                rec2=.false.
              endif
            enddo
            if(rec.and.vedg(rk2(n),1)>0)then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk2(n),1)
            endif
            if(rec2.and.vedg(rk2(n),2)>0.and.vedg(rk2(n),2)/=vedg(rk2(n),1))then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk2(n),2)
            endif
          else
            if(vedg(rk2(n),1)>0)then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk2(n),1)
            endif
            if(vedg(rk2(n),2)>0.and.vedg(rk2(n),2)/=vedg(rk2(n),1))then
              id=voronoiCell(cell)%vertexNb+1
              voronoiCell(cell)%vertexNb=id
              voronoiCell(cell)%vertexID(id)=vedg(rk2(n),2)
            endif
          endif
          if(vedg(rk2(n),1)<=0) voronoiCell(cell)%border=1

        else
          exit lp_ed2
        endif
      enddo lp_ed2

      ! Define the voronoi polygon as a convex hull
      if(voronoiCell(cell)%border==0)then
        vsort=-1
        n=voronoiCell(cell)%vertexNb
        do p=1,n
          id=voronoiCell(cell)%vertexID(p)
          vnx(p)=vcoordX(id)
          vny(p)=vcoordY(id)
        enddo

        ! Perform vertex sorting
        call Envelope(vnx(1:n),vny(1:n),n,vsort,nvert)

        ! Reorganise voronoi cell vertex
        if(nvert/=voronoiCell(cell)%vertexNb)then
          print*,nvert,voronoiCell(cell)%vertexNb
          print*,'Failed during creation of voronoi convex hull'
          call mpi_finalize(rc)
        endif

        do p=1,nvert
          sortedID(p)=voronoiCell(cell)%vertexID(vsort(p))
        enddo

        ! Update sorted vertex
        do p=1,nvert
          voronoiCell(cell)%vertexID(p)=sortedID(p)
        enddo

        ! Get voronoi cell parameters
        do p=1,nvert
          id=sortedID(p)
          id3=sortedID(p+1)
          if(p==nvert) id3=sortedID(1)
          ! Perimeter
          dist=sqrt((vcoordX(id3)-vcoordX(id))**2.0 &
              +(vcoordY(id3)-vcoordY(id))**2.0)
          voronoiCell(cell)%perimeter=voronoiCell(cell)%perimeter+dist
          ! Area
          area=(vcoordX(id)-tcoordX(cell))*(vcoordY(id3)-tcoordY(cell)) &
              -(vcoordX(id3)-tcoordX(cell))*(vcoordY(id)-tcoordY(cell))
          voronoiCell(cell)%area=voronoiCell(cell)%area+0.5_8*abs(area)
        enddo

      endif

    enddo

    ! Allocate voronoi edge length to delaunay neighborhood
    do cell=1,dnodes
       tsort=-1
       delaunayVertex(cell)%sortedHull=-1
       delaunayVertex(cell)%sortedHullNb=0
       if(voronoiCell(cell)%border==0)then
          vorPts=voronoiCell(cell)%vertexID

          do n=1,delaunayVertex(cell)%ngbNb
             id=delaunayVertex(cell)%ngbID(n)
             vorPtsn=voronoiCell(id)%vertexID
             m=0
             pt=-1
             do p=1,voronoiCell(cell)%vertexNb
                do k=1,voronoiCell(id)%vertexNb
                   if(vorPts(p)==vorPtsn(k))then
                      m=m+1
                      pt(m)=vorPts(p)
                   endif
                enddo
             enddo
             if(m==1)then
                dist=0.0_8
             elseif(m==0 .or. m>2)then
                print*,'Failed to allocate voronoi edge lenght to delaunay neighborhood'
                call mpi_finalize(rc)
             elseif(m==2)then
                dist=sqrt((vcoordX(pt(1))-vcoordX(pt(2)))**2.0+( &
                     vcoordY(pt(2))-vcoordY(pt(1)))**2.0)
             endif
             delaunayVertex(cell)%voronoi_edge(n)=dist
             delaunayVertex(cell)%distance(n)=sqrt((tcoordX(cell)-tcoordX(id))**2.0+( &
                tcoordY(cell)-tcoordY(id))**2.0)
             tnx(n)=tcoordX(id)
             tny(n)=tcoordY(id)
          enddo

          ! Perform vertex sorting
          n=delaunayVertex(cell)%ngbNb
          call Envelope(tnx(1:n),tny(1:n),n,tsort,tvert)

          delaunayVertex(cell)%sortedHullNb=tvert
          delaunayVertex(cell)%sortedHull=-1
          do p=1,tvert
            sortedID(p)=delaunayVertex(cell)%ngbID(tsort(p))
          enddo

          ! Update sorted vertex
          do p=1,tvert
            delaunayVertex(cell)%sortedHull(p)=sortedID(p)
          enddo

          if(tvert<delaunayVertex(cell)%ngbNb)then
            tt=0
            do while(tvert<delaunayVertex(cell)%ngbNb)
              sortedID=-1

              ! Find missing node ID
              loop:do p=1,delaunayVertex(cell)%ngbNb
                id=-1
                m=delaunayVertex(cell)%ngbID(p)
                do n=1,tvert
                  if(delaunayVertex(cell)%sortedHull(n)==m)id=0
                enddo
                if(id==-1)then
                  id=m
                  exit loop
                endif
              enddo loop

              if(id==0) goto 13
              if(tt>delaunayVertex(cell)%ngbNb*2) goto 13

              ! Add missing point in the convex hull
              lp:do p=1,tvert
                id1=delaunayVertex(cell)%sortedHull(p)
                if(p<tvert)then
                  id2=delaunayVertex(cell)%sortedHull(p+1)
                else
                  id2=delaunayVertex(cell)%sortedHull(1)
                endif
                if(PointTriangle(cell,id1,id2,id))then
                  do n=1,delaunayVertex(cell)%ngbNb
                    if(n<=p)then
                      sortedID(n)=delaunayVertex(cell)%sortedHull(n)
                      if(n==p)sortedID(n+1)=id
                    else
                      sortedID(n+1)=delaunayVertex(cell)%sortedHull(n)
                    endif
                  enddo
                  delaunayVertex(cell)%sortedHull=sortedID
                  tvert=tvert+1
                  exit lp
                endif
              enddo lp
              tt=tt+1
            enddo
13 continue
          endif
          do p=1,tvert
            delaunayVertex(cell)%sortedHull(p)=sortedID(p)
          enddo
        endif

      if(tcoordX(cell)==minx-dx.or.tcoordX(cell)==maxx+dx.or. &
        tcoordY(cell)==miny-dx.or.tcoordY(cell)==maxy+dx)then
        voronoiCell(cell)%border=1
      endif


      ! Count voronoi cells and vertices (used for output only)
      if(voronoiCell(cell)%border==0.and.voronoiCell(cell)%vertexNb>0)then
        vcellIN=vcellIN+1
        velemIN=velemIN+voronoiCell(cell)%vertexNb
      endif

    enddo

    return

  end subroutine DelaunayVoronoiDuality
  ! =====================================================================================
  function PointTriangle(id0,id1,id2,id) result(inside)

    logical::inside

    integer::id0,id1,id2,id
    real(kind=8)::denom,s,t,epsl,epsl2,xmin,xmax,ymin,ymax

    ! Point in triangle bounding box
    epsl=0.001
    epsl2=epsl*epsl
    xmin=min(tcoordX(id0),min(tcoordX(id1),tcoordX(id2)))-epsl
    xmax=max(tcoordX(id0),max(tcoordX(id1),tcoordX(id2)))+epsl
    ymin=min(tcoordY(id0),min(tcoordY(id1),tcoordY(id2)))-epsl
    ymax=max(tcoordY(id0),max(tcoordY(id1),tcoordY(id2)))+epsl

    if(tcoordX(id)<xmin.or.tcoordX(id)>xmax.or.tcoordY(id)<ymin.or.tcoordY(id)>ymax)then
      inside=.false.
    else
      inside=.false.
      denom=-tcoordY(id1)*tcoordX(id2)+tcoordY(id0)*(tcoordX(id2)-tcoordX(id1))+ &
        tcoordX(id0)*(tcoordY(id1)-tcoordY(id2))+tcoordX(id1)*tcoordY(id2)
      denom=1/denom
      s=tcoordY(id0)*tcoordX(id2)-tcoordX(id0)*tcoordY(id2)+(tcoordY(id2)-tcoordY(id0))*tcoordX(id)+ &
        (tcoordX(id0)-tcoordX(id2))*tcoordY(id)
      s=s*denom
      t=tcoordX(id0)*tcoordY(id1)-tcoordY(id0)*tcoordX(id1)+(tcoordY(id0)-tcoordY(id1))*tcoordX(id)+ &
        (tcoordX(id1)-tcoordX(id0))*tcoordY(id)
      t=t*denom
      if(s>0.and.t>0.0.and.1.0-s-t>0.0)then
        inside=.true.
      else
        if(distanceSquarePointToSegment(tcoordX(id0),tcoordY(id0), &
          tcoordX(id1),tcoordY(id1),tcoordX(id),tcoordY(id))<epsl2)then
          inside=.true.
        elseif(distanceSquarePointToSegment(tcoordX(id1),tcoordY(id1), &
          tcoordX(id2),tcoordY(id2),tcoordX(id),tcoordY(id))<epsl2)then
          inside=.true.
        elseif(distanceSquarePointToSegment(tcoordX(id2),tcoordY(id2), &
          tcoordX(id0),tcoordY(id0),tcoordX(id),tcoordY(id))<epsl2)then
          inside=.true.
        endif
      endif

    endif

  end function PointTriangle
  ! =====================================================================================
  function distanceSquarePointToSegment(x1,y1,x2,y2,x,y) result(dist2)

    real(kind=8)::x1,y1,x2,y2,x,y,sqrl,sqrl2,dotprod,dist2

    sqrl=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
    dotprod=((x -x1)*(x2-x1)+(y-y1)*(y2-y1))/sqrl

    if(dotprod<0.)then
      dist2=(x-x1)*(x-x1)+(y-y1)*(y-y1)
    elseif(dotprod<=1)then
      sqrl2=(x1-x)*(x1-x)+(y1-y)*(y1-y)
      dist2=sqrl2-dotprod*dotprod*sqrl
    else
      dist2=(x-x2)*(x-x2)+(y-y2)*(y-y2)
    endif

  end function distanceSquarePointToSegment
  ! =====================================================================================
  subroutine DelaunayBorders

    integer::k,p,n

    do k=1,dnodes
      voronoiCell(k)%btype=-1
      voronoiCell(k)%bpoint=-1
      if(voronoiCell(k)%border==1)then
        ! Corner points
        if(tcoordX(k)==minx-dx.and.tcoordY(k)==miny-dx)then
          voronoiCell(k)%btype=1 ! SW corner
          lp1: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==minx.and.tcoordY(n)==miny)then
              voronoiCell(k)%bpoint=n
              exit lp1
            endif
          enddo lp1
        elseif(tcoordX(k)==maxx+dx.and.tcoordY(k)==miny-dx)then
          voronoiCell(k)%btype=2 ! SE corner
          lp2: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==maxx.and.tcoordY(n)==miny)then
              voronoiCell(k)%bpoint=n
              exit lp2
            endif
          enddo lp2
        elseif(tcoordX(k)==minx-dx.and.tcoordY(k)==maxy+dx)then
          voronoiCell(k)%btype=3 ! NW corner
          lp3: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==minx.and.tcoordY(n)==maxy)then
              voronoiCell(k)%bpoint=n
              exit lp3
            endif
          enddo lp3
        elseif(tcoordX(k)==maxx+dx.and.tcoordY(k)==maxy+dx)then
          voronoiCell(k)%btype=4 ! NE corner
          lp4: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==maxx.and.tcoordY(n)==maxy)then
              voronoiCell(k)%bpoint=n
              exit lp4
            endif
          enddo lp4
        elseif(tcoordX(k)==minx-dx)then
          voronoiCell(k)%btype=5 ! West border
          lp5: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==minx.and.tcoordY(n)==tcoordY(k))then
              voronoiCell(k)%bpoint=n
              exit lp5
            endif
          enddo lp5
          if(voronoiCell(k)%bpoint<0)then
            lp5a: do p=1,delaunayVertex(k)%ngbNb
              n=delaunayVertex(k)%ngbID(p)
              if(tcoordX(n)==minx)then
                voronoiCell(k)%bpoint=n
                exit lp5a
              endif
            enddo lp5a
          endif
        elseif(tcoordX(k)==maxx+dx)then
          voronoiCell(k)%btype=6 ! East border
          lp6: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordX(n)==maxx.and.tcoordY(n)==tcoordY(k))then
              voronoiCell(k)%bpoint=n
              exit lp6
            endif
          enddo lp6
          if(voronoiCell(k)%bpoint<0)then
            lp6a: do p=1,delaunayVertex(k)%ngbNb
              n=delaunayVertex(k)%ngbID(p)
              if(tcoordX(n)==maxx)then
                voronoiCell(k)%bpoint=n
                exit lp6a
              endif
            enddo lp6a
          endif
        elseif(tcoordY(k)==miny-dx)then
          voronoiCell(k)%btype=7 ! South border
          lp7: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordY(n)==minx.and.tcoordX(n)==tcoordX(k))then
              voronoiCell(k)%bpoint=n
              exit lp7
            endif
          enddo lp7
          if(voronoiCell(k)%bpoint<0)then
            lp7a: do p=1,delaunayVertex(k)%ngbNb
              n=delaunayVertex(k)%ngbID(p)
              if(tcoordY(n)==miny)then
                voronoiCell(k)%bpoint=n
                exit lp7a
              endif
            enddo lp7a
          endif
        elseif(tcoordY(k)==maxy+dx)then
          voronoiCell(k)%btype=8 ! North border
          lp8: do p=1,delaunayVertex(k)%ngbNb
            n=delaunayVertex(k)%ngbID(p)
            if(tcoordY(n)==maxx.and.tcoordX(n)==tcoordX(k))then
              voronoiCell(k)%bpoint=n
              exit lp8
            endif
          enddo lp8
          if(voronoiCell(k)%bpoint<0)then
            lp8a: do p=1,delaunayVertex(k)%ngbNb
              n=delaunayVertex(k)%ngbID(p)
              if(tcoordY(n)==maxy)then
                voronoiCell(k)%bpoint=n
                exit lp8a
              endif
            enddo lp8a
          endif
        endif
        if(voronoiCell(k)%bpoint<0)then
          if(pet_id==0)then
            print*,'grid',minx-dx,maxx+dx,miny-dx,maxy+dx
            print*,'X Y',k,tcoordX(k),tcoordY(k)
            print*,'Border',k,voronoiCell(k)%btype,voronoiCell(k)%bpoint,voronoiCell(k)%border
            do p=1,delaunayVertex(k)%ngbNb
              n=delaunayVertex(k)%ngbID(p)
              print*,'ngh',p,tcoordX(n),tcoordY(n)
            enddo
          endif
        endif
      endif
    enddo

  end subroutine DelaunayBorders
  ! =====================================================================================
  subroutine Envelope(x,y,n,vertex,nvert)

    ! Find the vertices (in clockwise order) of a polygon enclosing
    ! the points (x(i), y(i), i=1, ..., n.
    ! On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
    ! iwk() is an integer work array which must have dimension at least n
    ! in the calling program.

    integer::n,vertex(n),nvert,iwk(n)
    integer::next(20),i,i1,i2,j,jp1,jp2,i2save,i3,i2next

    real(kind=8)::x(n),y(n),xmax,xmin,ymax,ymin,dist,dmax,dmin,x1,y1
    real(kind=8)::dx,dy,x2,y2,dx1,dx2,dmax1,dmax2,dy1,dy2,temp,zero

    zero=0.0_8

    if(n<2) return

    ! Choose the points with smallest & largest x- values as the
    ! first two vertices of the polygon.

    if(x(1)>x(n))then
       vertex(1)=n
       vertex(2)=1
       xmin=x(n)
       xmax=x(1 )
    else
       vertex(1)=1
       vertex(2)=n
       xmin=x(1)
       xmax=x(n)
    endif

    do i=2,n-1
       temp=x(i)
       if(temp<xmin)then
          vertex(1)=i
          xmin=temp
       elseif(temp>xmax)then
          vertex(2)=i
          xmax=temp
       endif
    enddo

    ! Special case, xmax = xmin
    if(xmax==xmin)then
       if(y(1)>y(n))then
          vertex(1)=n
          vertex(2)=1
          ymin=y(n)
          ymax=y(1)
       else
          vertex(1)=1
          vertex(2)=n
          ymin=y(1)
          ymax=y(n)
       endif

       do i=2,n-1
          temp=y(i)
          if(temp<ymin)then
             vertex(1)=i
             ymin=temp
          elseif(temp>ymax)then
             vertex(2)=i
             ymax=temp
          endif
       enddo

       nvert=2
       if(ymax==ymin) nvert=1
       return
    endif

    ! Set up two initial lists of points; those points above & those below the
    ! line joining the first two vertices. next(i) will hold the pointer to the
    ! point furthest from the line joining vertex(i) to vertex(i+1) on the left
    ! hand side.
    i1=vertex(1)
    i2=vertex(2)
    iwk(i1)=-1
    iwk(i2)=-1
    dx=xmax-xmin
    y1=y(i1)
    dy=y(i2)-y1
    dmax=zero
    dmin=zero
    next(1)=-1
    next(2)=-1

    do i=1,n
       if(i==vertex(1) .or. i==vertex(2)) cycle
       dist=(y(i)-y1)*dx-(x(i)-xmin)*dy
       if(dist>zero)then
          iwk(i1)=i
          i1=i
          if(dist>dmax)then
             next(1)=i
             dmax=dist
          endif
       elseif(dist<zero)then
          iwk(i2)=i
          i2=i
          if(dist<dmin)then
             next(2)=i
             dmin=dist
          endif
       endif
    enddo

    ! Ends of lists are indicated by pointers to -ve positions.
    iwk(i1)=-1
    iwk(i2)=-1
    nvert=2

    j=1

    ! Start of main process.
    ! Introduce new vertex between vertices j & j+1, if one has been found.
    ! Otherwise increase j. Exit if no more vertices.
40  if(next(j)<0)then
       if(j==nvert) return
       j=j+1
       goto 40
    endif

    jp1=j+1
    do i=nvert,jp1,-1
       vertex(i+1)=vertex(i)
       next(i+1)=next(i)
    enddo
    jp2=jp1+1
    nvert=nvert+1
    if(jp2>nvert) jp2=1
    i1=vertex(j)
    i2=next(j)
    i3=vertex(jp2)
    vertex(jp1)=i2

    ! Process the list of points associated with vertex j.   New list at vertex j
    ! consists of those points to the left of the line joining it to the new
    ! vertex (j+1).   Similarly for the list at the new vertex.
    ! Points on or to the right of these lines are dropped.
    x1=x(i1)
    x2=x(i2)
    y1=y(i1)
    y2=y(i2)
    dx1=x2-x1
    dx2=x(i3)-x2
    dy1=y2-y1
    dy2=y(i3)-y2
    DMAX1=zero
    dmax2=zero
    next(j)=-1
    next(jp1)=-1
    i2save=i2
    i2next=iwk(i2)
    i=iwk(i1)
    iwk(i1)=-1
    iwk(i2)=-1

60  if(i/=i2save)then
       dist=(y(i)-y1)*dx1-(x(i)-x1)*dy1
       if(dist>zero)then
          iwk(i1)=i
          i1=i
          if(dist>DMAX1)then
             next(j)=i
             DMAX1=dist
          endif
       else
          dist=(y(i)-y2)*dx2-(x(i)-x2)*dy2
          if(dist>zero)then
             iwk(i2)=i
             i2=i
             if(dist>dmax2)then
                next(jp1)=i
                dmax2=dist
             endif
          endif
       endif
       i=iwk(i)
    else
       i=i2next
    endif

    ! Get next point from old list at vertex j.
    if(i>0) goto 60

    ! End lists with -ve values.
    iwk(i1)=-1
    iwk(i2)=-1

    goto 40

  end subroutine Envelope
  ! =====================================================================================
  subroutine UnstructuredMeshDestroy

    if(allocated(dcentroid))deallocate(dcentroid)
    if(allocated(vedg))deallocate(vedg)
    if(allocated(delmt))deallocate(delmt)
    if(allocated(dedg))deallocate(dedg)
    if(allocated(tbound))deallocate(tbound)
    if(allocated(mask))deallocate(mask)
    if(allocated(unodeID))deallocate(unodeID)
    if(allocated(uelemID))deallocate(uelemID)
    if(allocated(delemID))deallocate(delemID)
    if(allocated(uownEID))deallocate(uownEID)
    if(allocated(outelem))deallocate(outelem)
    if(allocated(dglbID))deallocate(dglbID)
    if(allocated(outnode))deallocate(outnode)
    if(allocated(uownedID))deallocate(uownedID)
    if(allocated(unodeLID))deallocate(unodeLID)
    if(allocated(spmIDs))deallocate(spmIDs)
    if(allocated(elemtmask))deallocate(elemtmask)
    if(allocated(tcoordX))deallocate(tcoordX)
    if(allocated(tcoordY))deallocate(tcoordY)
    if(allocated(tcoordZ))deallocate(tcoordZ)
    if(allocated(tvertDisp))deallocate(tvertDisp)
    if(allocated(tflex))deallocate(tflex)
    if(allocated(gtflex))deallocate(gtflex)
    if(allocated(sedthick))deallocate(sedthick)
    if(allocated(ulay_th))deallocate(ulay_th)
    if(allocated(ulay_phi))deallocate(ulay_phi)
    if(allocated(ice_V))deallocate(ice_V)
    if(allocated(ice_H))deallocate(ice_H)
    if(allocated(vcoordX))deallocate(vcoordX)
    if(allocated(vcoordY))deallocate(vcoordY)
    if(allocated(refine_grid))deallocate(refine_grid)
    if(allocated(delaunayVertex))deallocate(delaunayVertex)
    if(allocated(voronoiCell))deallocate(voronoiCell)
    if(allocated(sedchange))deallocate(sedchange)
    if(allocated(lastsedthick))deallocate(lastsedthick)
    if(allocated(spmH))deallocate(spmH)
    if(allocated(Qs_in))deallocate(Qs_in)
    if(allocated(tempthick))deallocate(tempthick)

  end subroutine UnstructuredMeshDestroy
  ! =====================================================================================
  subroutine RegularGridDestroy

    if(allocated(bilinearX))deallocate(bilinearX)
    if(allocated(bilinearY))deallocate(bilinearY)
    if(allocated(bilinearV))deallocate(bilinearV)
    if(allocated(coordX))deallocate(coordX)
    if(allocated(coordY))deallocate(coordY)
    if(allocated(coordZ))deallocate(coordZ)
    if(allocated(snodeID))deallocate(snodeID)
    if(allocated(snodeLID))deallocate(snodeLID)
    if(allocated(oceanIDs))deallocate(oceanIDs)
    if(allocated(earthIDs))deallocate(earthIDs)
    if(allocated(selemID))deallocate(selemID)
    if(allocated(sownedID))deallocate(sownedID)
    if(allocated(rcoordX))deallocate(rcoordX)
    if(allocated(rcoordY))deallocate(rcoordY)
    if(allocated(rcoordZ))deallocate(rcoordZ)
    if(allocated(rtectoZ))deallocate(rtectoZ)
    if(allocated(rDisp))deallocate(rDisp)
    if(allocated(rhyDisp))deallocate(rhyDisp)
    if(allocated(rhxDisp))deallocate(rhxDisp)
    if(allocated(rvertDisp))deallocate(rvertDisp)
    if(allocated(rsedload))deallocate(rsedload)
    if(allocated(rainVal))deallocate(rainVal)
    if(allocated(frainmap))deallocate(frainmap)
    if(allocated(frainval))deallocate(frainval)
    if(allocated(rain_tstart))deallocate(rain_tstart)
    if(allocated(rain_tend))deallocate(rain_tend)
    if(allocated(fgeodyn))deallocate(fgeodyn)
    if(allocated(disp_time))deallocate(disp_time)
    if(allocated(vdisp_time))deallocate(vdisp_time)
    if(allocated(disp_fill))deallocate(disp_fill)
    if(allocated(vdisp_fill))deallocate(vdisp_fill)
    if(allocated(relmt))deallocate(relmt)
    if(allocated(sealvl))deallocate(sealvl)
    if(allocated(elaval))deallocate(elaval)

  end subroutine RegularGridDestroy
  ! =====================================================================================
end module geomesh
! =====================================================================================
