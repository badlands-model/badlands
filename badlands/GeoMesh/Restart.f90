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
!       Filename:  Restart.f90
!
!    Description:  Write SPM mesh HDF5
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module restart

  use hdf5
  use parameters
  use topology
  use hydroUtil
  use parallel
  use external_forces
  use kdtree2_module
  use kdtree2_precision_module

  implicit none


  integer::newnds,rstnodes,rstlay

  character(len=128)::frspm

  real(kind=8),dimension(:,:),allocatable::newcoord,rstXYZ,rstFphi,rstFth
  real(kind=8),dimension(:),allocatable::rstFload

contains

  ! =====================================================================================
  subroutine getSPM_hdf5topography3D

    logical::compression,simple

    integer::id,i,rank,rpet,localNd,p,id2,n,loclay

    real(kind=8),dimension(:),allocatable::nodes,nID,sedID,sedload,sedphi,sedth,pSedLoad
    real(kind=8),dimension(:,:),allocatable::prevNd,pFlexH,pFlexP
    real(kind=8)::x,y

    character(len=128)::text

    integer(hid_t)::file_id,d_spc
    integer(hid_t)::dset_id,dtype_id
    integer(hsize_t),dimension(2)::dims,maxdims

    rstnodes=0

    if(allocated(ulay_th)) deallocate(ulay_th)
    if(allocated(ulay_phi)) deallocate(ulay_phi)

    do rpet=0,restartPet-1
        frspm=rstfolder
        call noblnk(frspm)
        text='/outputs/SurfaceSPM.'
        call append_str(frspm,text)
        call append_zero(frspm,restartStep)
        text='.p'
        call append_str(frspm,text)
        call append_nb(frspm,rpet)
        text='.h5'
        call append_str(frspm,text)

        ! Initialize predefined datatypes
        call h5open_f(rc)
        call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

        ! Open the file collectively.
        call h5fopen_f(frspm,h5f_acc_rdonly_f,file_id,rc)

        ! The node global ID
        text="/nID"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nID(dims(1)*dims(2)))
        localNd=dims(2)
        call h5dread_f(dset_id,h5t_native_double,nID,dims,rc)

        ! The Coordinates - vertices
        text="/vertices"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nodes(3*localNd))
        if(rpet==0)then
            allocate(prevNd((restartPet+1)*localNd,4))
        endif
        call h5dread_f(dset_id,h5t_native_double,nodes,dims,rc)

        ! The sediment pile thickness
        text="/cumdz"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(sedID(dims(1)*dims(2)))
        call h5dread_f(dset_id,h5t_native_double,sedID,dims,rc)

        ! Flexure data
        if(flexure)then
          text="/flex"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedload(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedloader,dims,rc)

          text="/layth"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedth(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedth,dims,rc)

          text="/layphi"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedphi(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedphi,dims,rc)
          ! Update number of stratigrahic layers
          loclay=dims(1)
          if(rpet==0)then
            allocate(pSedLoad((restartPet+1)*localNd))
            allocate(pFlexH((restartPet+1)*localNd,loclay))
            allocate(pFlexP((restartPet+1)*localNd,loclay))
          endif
        endif

        id=1
        id2=0
        do p=1,localNd
            i=int(nID(p))
            rstnodes=max(rstnodes,i)
            prevNd(i,1)=nodes(id)
            prevNd(i,2)=nodes(id+1)
            prevNd(i,3)=nodes(id+2)
            prevNd(i,4)=sedID(p)
            id=id+3
            if(flexure)then
              pSedLoad(i)=sedload(p)
              do n=1,loclay
                id2=id2+1
                pFlexH(i,n)=sedth(id2)
                pFlexP(i,n)=sedphi(id2)
              enddo
            endif
        enddo

        ! Close the dataset
        call h5sclose_f(d_spc,rc)
        call h5dclose_f(dset_id,rc)
        ! Close the file.
        call h5fclose_f(file_id,rc)
        ! Close interface
        call h5close_f(rc)
        deallocate(nodes,nID,sedID,sedphi,sedth,sedload)
    enddo

    if(allocated(rstXYZ)) deallocate(rstXYZ)
    if(allocated(newcoord)) deallocate(newcoord)
    allocate(rstXYZ(rstnodes,4))
    allocate(newcoord(rstnodes,2))
    if(flexure)then
      rstlay=loclay
      if(allocated(rstFphi)) deallocate(rstFphi)
      if(allocated(rstFth)) deallocate(rstFth)
      if(allocated(rstFload)) deallocate(rstFload)
      allocate(rstFphi(rstnodes,loclay))
      allocate(rstFth(rstnodes,loclay))
      allocate(rstFload(rstnodes))
    endif

    newnds=0
    do id=1,rstnodes
        x=prevNd(id,1)
        y=prevNd(id,2)
        rstXYZ(id,1:4)=prevNd(id,1:4)
        if(x>minx+dx.and.x<maxx-dx.and.y>miny+dx.and.y<maxy-dx)then
            newnds=newnds+1
            newcoord(newnds,1)=x
            newcoord(newnds,2)=y
        endif
        if(flexure)then
          rstFload(id)=pSedLoad(id)
          do n=1,rstlay
            rstFth(id,n)=pFlexH(id,n)
            rstFphi(id,n)=pFlexP(id,n)
          enddo
        endif
    enddo

    deallocate(prevNd,pSedLoad,pFlexH,pFlexP)

    return

  end subroutine getSPM_hdf5topography3D
  ! =====================================================================================
  subroutine getRestart_topography

    integer::id,k,rc,p,n,step

    real(kind=8),dimension(2)::txy
    real(kind=8),dimension(2,rstnodes)::Fd
    real(kind=8),dimension(:),allocatable::Fphi,Fth

    type(kdtree2),pointer::Ftree
    type(kdtree2_result),dimension(1)::FRslt

    if(flexure)then
      if(allocated(sedloader)) deallocate(sedloader)
      if(allocated(ulay_phi)) deallocate(ulay_phi)
      if(allocated(ulay_th)) deallocate(ulay_th)
      if(flex_dx<dx) flex_dx=dx
      step=int(flex_dx/dx)
      nbfx=int(nx/step)
      nbfy=int(ny/step)
      if(allocated(prevload)) deallocate(prevload)
      allocate(prevload(nbfx+2,nbfy+2))
      if(allocated(sedloader)) deallocate(sedloader)
      allocate(sedloader((nbfx+2)*(nbfy+2)))
      flex_lay=flex_lay+rstlay
      allocate(ulay_th(dnodes,flex_lay))
      allocate(ulay_phi(dnodes,flex_lay))
      allocate(Fphi(dnodes*rstlay))
      allocate(Fth(dnodes*rstlay))
    endif

    if(pet_id==0)then
        do k=1,rstnodes
            Fd(1,k)=rstXYZ(k,1)
            Fd(2,k)=rstXYZ(k,2)
        enddo
        Ftree=>kdtree2_create(Fd,sort=.true.,rearrange=.true.)

        n=0
        do k=1,dnodes
            txy(1)=tcoordX(k)
            txy(2)=tcoordY(k)
            call kdtree2_n_nearest(Ftree,txy,nn=1,results=FRslt)
            tcoordZ(k)=rstXYZ(FRslt(1)%idx,3)
            sedthick(k)=rstXYZ(FRslt(1)%idx,4)
            if(tcoordX(k)>minx.and.tcoordX(k)<maxx &
              .and.tcoordY(k)>miny.and.tcoordY(k)<maxy)then
              sedthick(k)=sedthick(k)+100000.
            else
              sedthick(k)=100000.
            endif
            if(flexure)then
              sedloader(k)=rstFload(FRslt(1)%idx)
              do p=1,rstlay
                n=n+1
                Fth(n)=rstFth(FRslt(1)%idx,p)
                Fphi(n)=rstFphi(FRslt(1)%idx,p)
              enddo
            endif
        enddo
        call kdtree2_destroy(Ftree)
    endif

    if(allocated(bilinearX)) deallocate(bilinearX)
    if(allocated(bilinearY)) deallocate(bilinearY)
    if(allocated(bilinearV)) deallocate(bilinearV)
    allocate(bilinearX(nx+2))
    allocate(bilinearY(ny+2))
    allocate(bilinearV(nx+2,ny+2))
    if(disp3d)then
        if(allocated(bilinearHx)) deallocate(bilinearHx)
        if(allocated(bilinearHy)) deallocate(bilinearHy)
        allocate(bilinearHx(nx+2,ny+2))
        allocate(bilinearHy(nx+2,ny+2))
    endif
    bilinearX(1)=real(minx-dx)
    do id=2,nx+2
        bilinearX(id)=bilinearX(id-1)+real(dx)
    enddo
    bilinearY(1)=real(miny-dx)
    do id=2,ny+2
        bilinearY(id)=bilinearY(id-1)+real(dx)
    enddo

    if(allocated(rstXYZ)) deallocate(rstXYZ)
    if(flexure)then
      if(allocated(rstFth)) deallocate(rstFth)
      if(allocated(rstFphi)) deallocate(rstFphi)
      if(allocated(rstFload)) deallocate(rstFload)
    endif
    call mpi_bcast(tcoordZ,dnodes,mpi_double_precision,0,badlands_world,rc)
    call mpi_bcast(sedthick,dnodes,mpi_double_precision,0,badlands_world,rc)

    if(flexure)then
      call mpi_bcast(sedloader,(nbfx+2)*(nbfy+2),mpi_double_precision,0,badlands_world,rc)
      call mpi_bcast(Fth,dnodes*rstlay,mpi_double_precision,0,badlands_world,rc)
      call mpi_bcast(Fphi,dnodes*rstlay,mpi_double_precision,0,badlands_world,rc)
      flex_lay=rstlay
      n=0
      do k=1,dnodes
        do p=1,rstlay
          n=n+1
          ulay_th(k,p)=Fth(n)
          ulay_phi(k,p)=Fphi(n)
        enddo
      enddo
      deallocate(Fth,Fphi)
    endif

    return

  end subroutine getRestart_topography
  ! =====================================================================================
  subroutine getSPM_hdf5topography

    logical::compression,simple

    integer::id,i,rank,rpet,localNd,p,n,loclay,id2,step,j

    real(kind=8),dimension(:),allocatable::nodes,nID,sedID,sedlID,layphi,layth,sedload, cflex

    character(len=128)::text

    integer(hid_t)::file_id,d_spc
    integer(hid_t)::dset_id,dtype_id
    integer(hsize_t),dimension(2)::dims,maxdims

    if(flexure)then
      if(allocated(sedloader)) deallocate(sedloader)
      if(allocated(ulay_phi)) deallocate(ulay_phi)
      if(allocated(ulay_th)) deallocate(ulay_th)
      !allocate(sedloader(dnodes))

      if(flex_dx<dx) flex_dx=dx
      step=int(flex_dx/dx)
      nbfx=int(nx/step)
      nbfy=int(ny/step)
      if(allocated(prevload)) deallocate(prevload)
      allocate(prevload(nbfx+2,nbfy+2))
      if(allocated(sedloader)) deallocate(sedloader)
      allocate(sedloader((nbfx+2)*(nbfy+2)))
    endif
    !if(allocated(ulay_th)) deallocate(ulay_th)
    !if(allocated(ulay_phi)) deallocate(ulay_phi)

    do rpet=0,restartPet-1
        frspm=rstfolder
        call noblnk(frspm)
        text='/outputs/SurfaceSPM.'
        call append_str(frspm,text)
        call append_zero(frspm,restartStep)
        text='.p'
        call append_str(frspm,text)
        call append_nb(frspm,rpet)
        text='.h5'
        call append_str(frspm,text)

        ! Initialize predefined datatypes
        call h5open_f(rc)
        call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

        ! Open the file collectively.
        call h5fopen_f(frspm,h5f_acc_rdonly_f,file_id,rc)

        ! The node global ID
        text="/nID"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nID(dims(1)*dims(2)))
        localNd=dims(2)
        call h5dread_f(dset_id,h5t_native_double,nID,dims,rc)

        ! The Coordinates - vertices
        text="/vertices"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nodes(3*localNd))

        call h5dread_f(dset_id,h5t_native_double,nodes,dims,rc)

        ! The sediment pile thickness
        text="/cumdz"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(sedID(dims(1)*dims(2)))
        call h5dread_f(dset_id,h5t_native_double,sedID,dims,rc)

        ! The last sediment pile thickness for flexure
        if(flexure)then
          text="/flex"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedload(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedload,dims,rc)

          text="/cumflex"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(cflex(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,cflex,dims,rc)

          text="/layth"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(layth(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,layth,dims,rc)

          text="/layphi"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(layphi(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,layphi,dims,rc)
          ! Update number of stratigrahic layers
          loclay=dims(1)
          flex_lay=flex_lay+dims(1)
          if(.not.allocated(ulay_th))then
            allocate(ulay_th(dnodes,flex_lay),ulay_phi(dnodes,flex_lay))
            !ulay_th=0.
            !ulay_phi=0.
          endif
          if(.not.allocated(gtflex))then
            allocate(gtflex(dnodes)) 
            !gtflex = 0.
          endif
        endif

        if(flexure)then
          p = 0
          do j=1,nbfy+2
            do i=1,nbfx+2
              p=p+1
              prevload(i,j)=sedload(p)
              sedloader(p)=sedload(p)
            enddo
          enddo
        endif

        id=1
        id2=0
        do p=1,localNd
            i=int(nID(p))
            if(abs(nodes(id)-tcoordX(i))>0.01.or.abs(nodes(id+1)-tcoordY(i))>0.01)then
                print*,'Problem during restart',i,nodes(id),tcoordX(i),nodes(id+1),tcoordY(i)
                call mpi_finalize(rc)
            endif
            tcoordZ(i)=nodes(id+2)
            if(tcoordX(i)>minx.and.tcoordX(i)<maxx &
              .and.tcoordY(i)>miny.and.tcoordY(i)<maxy)then
              sedthick(i)=sedID(p)+100000.
            else
              sedthick(i)=100000.
            endif
            id=id+3
            if(flexure)then
              gtflex(i) = cflex(p)
              do n=1,loclay
                id2=id2+1
                ulay_th(i,n)=layth(id2)
                ulay_phi(i,n)=layphi(id2)
              enddo
            endif
        enddo

        ! Close the dataset
        call h5sclose_f(d_spc,rc)
        call h5dclose_f(dset_id,rc)

        ! Close the file.
        call h5fclose_f(file_id,rc)

        ! Close interface
        call h5close_f(rc)

        deallocate(nodes,nID,sedID)
        if( allocated(sedlid) ) deallocate(sedlID)
        if(flexure)deallocate(layth,layphi,sedload,cflex)
    enddo

    if(allocated(bilinearX)) deallocate(bilinearX)
    if(allocated(bilinearY)) deallocate(bilinearY)
    if(allocated(bilinearV)) deallocate(bilinearV)
    allocate(bilinearX(nx+2))
    allocate(bilinearY(ny+2))
    allocate(bilinearV(nx+2,ny+2))
    if(disp3d)then
        if(allocated(bilinearHx)) deallocate(bilinearHx)
        if(allocated(bilinearHy)) deallocate(bilinearHy)
        allocate(bilinearHx(nx+2,ny+2))
        allocate(bilinearHy(nx+2,ny+2))
    endif
    bilinearX(1)=real(minx-dx)
    do id=2,nx+2
        bilinearX(id)=bilinearX(id-1)+real(dx)
    enddo
    bilinearY(1)=real(miny-dx)
    do id=2,ny+2
        bilinearY(id)=bilinearY(id-1)+real(dx)
    enddo

    flex_lay=loclay

    return

  end subroutine getSPM_hdf5topography
  ! =====================================================================================
end module restart
