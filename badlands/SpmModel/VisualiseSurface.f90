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
!       Filename:  VisualiseSurface.f90
!
!    Description:  Implements the XdmF and HdF5 SPM output of SPM surface evolution
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module outspm_surface

  use hdf5
  use parallel
  use parameters
  use topology
  use FoX_wxml
  use hydroUtil
  use hydrology
  use external_forces

  implicit none

  character(len=128)::fspm

contains

  ! =====================================================================================
  subroutine spm_hdf5(iter)

    logical::compression

    integer::id,id2,i,j,k,p,rank,iter,totnodes,totelems,ierr
    integer,dimension(:),allocatable::connect
    real(kind=8),dimension(:),allocatable::nodes,facc,dz,cumdz,rego,cID,nID,sl
    real(kind=8),dimension(:),allocatable::coldU,coldH,sedflex,cumflex
    real(kind=8),dimension(:),allocatable::sedth,sedphi

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    fspm='SurfaceSPM'

    file=''
    file=fspm
    call noblnk(file)
    text='.'
    call append_str(file,text)
    call append_zero(file,iter)
    text='.p'
    call append_str(file,text)
    call append_nb(file,pet_id)
    text='.h5'
    call append_str(file,text)
    call addpath1(file)
    totnodes=upartN
    totelems=delemoo

    allocate(nodes(3*totnodes))
    allocate(dz(totnodes))
    allocate(cumdz(totnodes))
    allocate(rego(totnodes))
    allocate(facc(totnodes))
    allocate(cID(totnodes))
    allocate(nID(totnodes))
    allocate(sl(totnodes))
    allocate(connect(3*totelems))
    if(ice_dx>0.)then
      allocate(coldU(totnodes))
      allocate(coldH(totnodes))
    endif
    if(flexure)then
      allocate(cumflex(totnodes))
      allocate(sedth(totnodes*flex_lay))
      allocate(sedphi(totnodes*flex_lay))
    endif

    ! Create nodes arrays
    id=1
    id2=0
    do k=1,totnodes
       i=unodeID(k)
       nodes(id)=tcoordX(i)
       nodes(id+1)=tcoordY(i)
       if(voronoiCell(i)%border==1)then
         nodes(id+2)=spmZ(voronoiCell(i)%bpoint)
       else
         nodes(id+2)=spmZ(i)
       endif
       if(voronoiCell(i)%border==1)then
        facc(k)=0.
        dz(k)=0.0
        rego(k)=0.0
       else
        facc(k)=discharge(i)
        if(watercell(i)<0.)then
            rego(k)=0.
        else
            rego(k)=watercell(i)
        endif
       endif
       if(tcoordX(i)>minx.and.tcoordX(i)<maxx &
          .and.tcoordY(i)>miny.and.tcoordY(i)<maxy)then
          cumdz(k)=sedthick(i)-100000.
          dz(k)=sedthick(i)-lastsedthick(i)
       else
          cumdz(k)=0.
          dz(k)=0.0
       endif
       if(flexure)then
         cumflex(k)=gtflex(i)
         do p=1,flex_lay
           id2=id2+1
           sedth(id2)=ulay_th(i,p)
           sedphi(id2)=ulay_phi(i,p)
         enddo
       endif
       cID(k)=real(bsID(i))
       sl(k)=gsea%actual_sea
       if(ice_dx>0.)then
        if(tcoordX(i)>minx.and.tcoordX(i)<maxx &
            .and.tcoordY(i)>miny.and.tcoordY(i)<maxy)then
          coldU(k)=ice_V(i)
          coldH(k)=ice_H(i)
        else
          coldU(k)=0.0
          coldH(k)=0.0
        endif
       endif
       nID(k)=real(i)
       id=id+3
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f(h5p_file_access_f,plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    ! Surface - connectivity
    id=1
    p=0
    do k=1,totelems
        i=delemID(k)
        if(elemtmask(i)==0)then
          connect(id)=unodeLID(delmt(i,1))
          connect(id+1)=unodeLID(delmt(i,2))
          connect(id+2)=unodeLID(delmt(i,3))
          id=id+3
        endif
    enddo

    dims(1)=3
    dims(2)=totelems
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/connectivity"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    dims(1)=1
    dims(2)=totelems*3

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,connect,dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! Nodes ID
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/nID"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,nID,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! The Coordinates - vertices
    dims(1)=3
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/vertices"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc )
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    dims(1)=1
    dims(2)=totnodes*3

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,nodes,dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! Flow accumulation
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/facc"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,facc,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Cumulative change
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/cumdz"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,cumdz,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Elevation change
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/dz"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,dz,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Flexural thickness
    if(flexure)then
      ! Cumulative flexure
      dims(1)=1
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/cumflex"

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,cumflex,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

      ! Sediment load for flexural isostasy computation
      dims(1)=1
      dims(2)=(nbfx+2)*(nbfy+2)
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/flex"
      allocate(sedflex(dims(2)))
      p=0
      do j=1,nbfy+2
        do i=1,nbfx+2
          p=p+1
          sedflex(p)=prevload(i,j)
        enddo
      enddo

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,sedflex,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

      ! Sediment layer thickness
      dims(1)=flex_lay
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="layth"

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,sedth,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

      ! Sediment layer porosity
      dims(1)=flex_lay
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/layphi"

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,sedphi,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

    endif

    ! Catchment ID
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/cID"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,cID,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Sea-level
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/sl"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,sl,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    if(ice_dx>0.)then
      ! Ice thickness
      dims(1)=1
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/iceH"

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,coldH,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

      ! Ice flow
      dims(1)=1
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/iceU"

      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)

      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,coldU,dims,ierr)
      call h5pclose_f(plist_id,ierr)

      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)
    endif

    ! Regolith depth
    if(regoProd>0.)then
        dims(1)=1
        dims(2)=totnodes
        rank=2
        call h5screate_simple_f(rank,dims,filespace,ierr)
        text=''
        text="/reg"

        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
        call h5pset_deflate_f(plist_id,9,ierr)
        call h5pset_chunk_f(plist_id,rank,dims,ierr)

        ! Create the dataset with default properties
        call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

        ! Write the dataset collectively
        call h5dwrite_f(dset_id,h5t_native_double,rego,dims,ierr)
        call h5pclose_f(plist_id,ierr)

        ! Close the dataset
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(filespace,ierr)
    endif

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(nodes,connect,facc,dz,rego,cID,nID)
    if(allocated(cumflex))deallocate(cumflex)

    return

  end subroutine spm_hdf5
  ! =====================================================================================
  subroutine spm_xmf(iter)

    type(xmlf_t)::xf
    integer::ierr
    integer::iter,totnodes,totelems,k
    character(len=128)::str,stg,filename,filename1,filename2,file,filename3,filename10
    character(len=128)::filename4,filename5,filename6,filename7,filename8,filename9

    call spm_hdf5(iter)

    if(.not.allocated(outnode))then
        allocate(outnode(npets),outelem(npets))
        call mpi_gather(upartN,1,mpi_integer,outnode,1,mpi_integer,0,badlands_world,ierr)
        call mpi_gather(delemoo,1,mpi_integer,outelem,1,mpi_integer,0,badlands_world,ierr)
    endif

    if(pet_id==0)then
        fspm='SurfaceSPM'
        file=''
        file=fspm
        call noblnk(file)
        str='.'
        call append_str(file,str)
        call append_zero(file,iter)
        str='.xmf'
        call append_str(file,str)
        call addpath1(file)

        call xml_OpenFile(file,xf)

        ! Header
        call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
        call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
        call xml_NewElement(xf,"Xdmf")
        call xml_AddAttribute(xf,"Version","2.0")
        call xml_NewElement(xf,"Domain")
        call xml_NewElement(xf,"Grid")
        call xml_AddAttribute(xf,"GridType","Collection")
        call xml_AddAttribute(xf,"CollectionType","Spatial")
        call xml_NewElement(xf,"Time")
        call xml_AddAttribute(xf,"Type","Single")
        call xml_AddAttribute(xf,"Value",time_display)
        call xml_EndElement(xf,"Time")
        do k=1,npets
            totnodes=outnode(k)
            totelems=outelem(k)
            filename=''
            filename=fspm
            call noblnk(filename)
            str='.'
            call append_str(filename,str)
            call append_zero(filename,iter)
            str='.p'
            call append_str(filename,str)
            call append_nb(filename,k-1)
            str='.h5'
            call append_str(filename,str)
            filename1=filename
            filename2=filename
            filename3=filename
            filename4=filename
            filename5=filename
            filename6=filename
            filename7=filename
            filename8=filename
            filename9=filename
            filename10=filename
            str=':/connectivity'
            call append_str(filename,str)
            str=':/vertices'
            call append_str(filename1,str)
            str=':/facc'
            call append_str(filename2,str)
            str=':/dz'
            call append_str(filename3,str)
            str=':/reg'
            call append_str(filename4,str)
            str=':/cID'
            call append_str(filename5,str)
            str=':/sl'
            call append_str(filename6,str)
            str=':/iceH'
            call append_str(filename7,str)
            str=':/iceU'
            call append_str(filename8,str)
            str=':/cumdz'
            call append_str(filename9,str)
            str=':/cumflex'
            call append_str(filename10,str)

            ! Block begin
            call xml_NewElement(xf,"Grid")
            str='SurfBlock.'
            call append_zero(str,iter)
            stg='.p'
            call append_str(str,stg)
            call append_zero(str,k-1)
            call xml_AddAttribute(xf,"Name",trim(str))
            call xml_NewElement(xf,"Topology")
            call xml_AddAttribute(xf,"Type","Triangle")
            call xml_AddAttribute(xf,"NumberOfElements",totelems)
            call xml_AddAttribute(xf,"BaseOffset","1")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"DataType","Int")
            str=' '
            call append_nb2(str,totelems)
            call append_nb2(str,3)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Topology")

            ! Geometry
            call xml_NewElement(xf,"Geometry")
            call xml_AddAttribute(xf,"Type","XYZ")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","8")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,3)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename1))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Geometry")

            ! Flow accumulation
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Flow discharge")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,1)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename2))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Cumulative elevation change
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Cumulative elevation change")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,1)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename9))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            if(flexure)then
              !Cumulative flexural isostasy
              call xml_NewElement(xf,"Attribute")
              call xml_AddAttribute(xf,"Type","Scalar")
              call xml_AddAttribute(xf,"Center","Node")
              call xml_AddAttribute(xf,"Name","Cumulative flexural isostasy")
              call xml_NewElement(xf,"DataItem")
              call xml_AddAttribute(xf,"Format","HDF")
              call xml_AddAttribute(xf,"NumberType","Float")
              call xml_AddAttribute(xf,"Precision","4")
              str=' '
              call append_nb2(str,totnodes)
              call append_nb2(str,1)
              call xml_AddAttribute(xf,"Dimensions",trim(str))
              call xml_AddCharacters(xf,trim(filename10))
              call xml_EndElement(xf,"DataItem")
              call xml_EndElement(xf,"Attribute")
            endif

            ! Elevation change
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Elevation change")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,1)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename3))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Catchment ID
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Catchment ID")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,1)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename5))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Sea level
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Sealevel")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            str=' '
            call append_nb2(str,totnodes)
            call append_nb2(str,1)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename6))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            if(ice_dx>0.)then
              ! Ice thickness
              call xml_NewElement(xf,"Attribute")
              call xml_AddAttribute(xf,"Type","Scalar")
              call xml_AddAttribute(xf,"Center","Node")
              call xml_AddAttribute(xf,"Name","Ice thickness")
              call xml_NewElement(xf,"DataItem")
              call xml_AddAttribute(xf,"Format","HDF")
              call xml_AddAttribute(xf,"NumberType","Float")
              call xml_AddAttribute(xf,"Precision","4")
              str=' '
              call append_nb2(str,totnodes)
              call append_nb2(str,1)
              call xml_AddAttribute(xf,"Dimensions",trim(str))
              call xml_AddCharacters(xf,trim(filename7))
              call xml_EndElement(xf,"DataItem")
              call xml_EndElement(xf,"Attribute")

              ! Ice flow velocity
              call xml_NewElement(xf,"Attribute")
              call xml_AddAttribute(xf,"Type","Scalar")
              call xml_AddAttribute(xf,"Center","Node")
              call xml_AddAttribute(xf,"Name","Ice flow velocity")
              call xml_NewElement(xf,"DataItem")
              call xml_AddAttribute(xf,"Format","HDF")
              call xml_AddAttribute(xf,"NumberType","Float")
              call xml_AddAttribute(xf,"Precision","4")
              str=' '
              call append_nb2(str,totnodes)
              call append_nb2(str,1)
              call xml_AddAttribute(xf,"Dimensions",trim(str))
              call xml_AddCharacters(xf,trim(filename8))
              call xml_EndElement(xf,"DataItem")
              call xml_EndElement(xf,"Attribute")
            endif

            ! Regolith depth
            if(regoProd>0.)then
                call xml_NewElement(xf,"Attribute")
                call xml_AddAttribute(xf,"Type","Scalar")
                call xml_AddAttribute(xf,"Center","Node")
                call xml_AddAttribute(xf,"Name","Regolith depth")
                call xml_NewElement(xf,"DataItem")
                call xml_AddAttribute(xf,"Format","HDF")
                call xml_AddAttribute(xf,"NumberType","Float")
                call xml_AddAttribute(xf,"Precision","4")
                str=' '
                call append_nb2(str,totnodes)
                call append_nb2(str,1)
                call xml_AddAttribute(xf,"Dimensions",trim(str))
                call xml_AddCharacters(xf,trim(filename4))
                call xml_EndElement(xf,"DataItem")
                call xml_EndElement(xf,"Attribute")
            endif
            call xml_EndElement(xf,"Grid")
        enddo

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    return

  end subroutine spm_xmf
  ! =====================================================================================
  subroutine visualise_surface_changes(iter)

    ! Parameters Declaration
    type(xmlf_t)::xf

    integer::i,iter,it0
    character(len=128)::filename,str,fname

    call spm_xmf(iter)

    if(pet_id==0)then
        filename='SPMsurface_series.xdmf'
        call addpath1(filename)

        ! Header
        call xml_OpenFile(filename,xf)
        call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
        call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
        call xml_NewElement(xf,"Xdmf")
        call xml_AddAttribute(xf,"Version","2.0")
        call xml_NewElement(xf,"Domain")
        call xml_NewElement(xf,"Grid")
        call xml_AddAttribute(xf,"GridType","Collection")
        call xml_AddAttribute(xf,"CollectionType","Temporal")

        it0=1
        ! Loop over time step
        do i=it0,iter+1
            ! Grid name
            fname=''
            fname=fspm
            call noblnk(fname)
            str='.'
            call append_str(fname,str)
            call append_zero(fname,i-1)
            str='.xmf'
            call append_str(fname,str)
            call xml_NewElement(xf,"xi:include")
            call xml_AddAttribute(xf,"href",trim(fname))
            call xml_AddAttribute(xf,"xpointer","xpointer(//Xdmf/Domain/Grid)")
            call xml_EndElement(xf,"xi:include")
        enddo

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    return

   end subroutine visualise_surface_changes
   ! =====================================================================================
end module outspm_surface
