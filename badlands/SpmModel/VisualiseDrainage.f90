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
!       Filename:  VisualiseDrainage.f90
!
!    Description:  Implements the XdmF and HdF5 SPM output of SPM drainage evolution
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module outspm_drainage

  use hdf5
  use parameters
  use topology
  use FoX_wxml
  use hydroUtil
  use hydrology

  implicit none

  character(len=128)::fdspm

contains

  ! =====================================================================================
  subroutine drainage_hdf5(iter,totelems)

    logical::compression

    integer::id,i,k,p,lid,rank,iter,totnodes,totelems,ierr,rcv
    integer,dimension(:),allocatable::connect,drainers,connectg

    real(kind=8),dimension(:),allocatable::nodes,facc,cID,stNb,chiID

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    fdspm='DrainageSPM'

    file=''
    file=fdspm
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
    totnodes=drainOde !dnodes

    allocate(nodes(3*totnodes))
    allocate(facc(totnodes))
    allocate(cID(totnodes))
    allocate(chiID(totnodes))
    allocate(stNb(totnodes))
    allocate(drainers(totnodes))
    allocate(connect(2*totnodes))
    allocate(connectg(2*totnodes))
    if(.not.allocated(dglbID)) allocate(dglbID(dnodes))

    ! Create nodes arrays
    id=1
    p=0
    dglbID=-1
    drainers=-1
    facc=0.
    cID=0.
    stNb=0.
    do lid=1,localNodes !dnodes
       k=localNodesGID(lid)
       i=stackOrder(k)
       p=p+1
       nodes(id)=tcoordX(i)
       nodes(id+1)=tcoordY(i)
       if(voronoiCell(i)%border==1)then
         nodes(id+2)=spmZ(voronoiCell(i)%bpoint)
       else
         nodes(id+2)=spmZ(i)
       endif
       facc(p)=discharge(i)
       chiID(p)=chi(i)
       cID(p)=real(bsID(i))
       stNb(p)=real(strahler(i))
       id=id+3
       drainers(p)=i
       dglbID(i)=p
       if(subcatchmentProc(i)/= &
         subcatchmentProc(receivers(i)))then
         i=receivers(i)
         p=p+1
         nodes(id)=tcoordX(i)
         nodes(id+1)=tcoordY(i)
         if(voronoiCell(i)%border==1)then
            nodes(id+2)=spmZ(voronoiCell(i)%bpoint)
         else
            nodes(id+2)=spmZ(i)
         endif
         facc(p)=discharge(i)
         chiID(p)=chi(i)
         cID(p)=real(bsID(i))
         stNb(p)=real(strahler(i))
         id=id+3
         drainers(p)=i
         dglbID(i)=p
       endif
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f( h5p_file_access_f, plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    ! Drainage - connectivity
    id=1
    totelems=0
    p=0
    do k=1,drainOde
        i=drainers(k)
        p=p+1
        if(i>0)then
            if(voronoiCell(i)%border==0.and.voronoiCell(i)%btype<0.and.tcoordX(i)>=minx.and.&
             tcoordX(i)<=maxx.and.tcoordY(i)>=miny.and.tcoordY(i)<=maxy)then
                rcv=receivers(i)
                if(rcv/=i.and.voronoiCell(rcv)%border==0.and.voronoiCell(rcv)%btype<0)then
                    if(subcatchmentProc(i)==pet_id.and.tcoordX(rcv)>=minx.and.tcoordX(rcv)<=maxx.and. &
                      tcoordY(rcv)>=miny.and.tcoordY(rcv)<=maxy)then
                        totelems=totelems+1
                        if(totnodes*2<totelems)print*,'Problem when writing drainage hdf5 files'
                        connect(id)=k
                        connect(id+1)=dglbID(rcv)
                        connectg(id)=i
                        connectg(id+1)=rcv
                        id=id+2
                    endif
                endif
            endif
        endif
    enddo

    dims(1)=2
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
    dims(2)=totelems*2

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,connect(1:totelems*2),dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! Global connectivity
    dims(1)=2
    dims(2)=totelems
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/globconnect"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    dims(1)=1
    dims(2)=totelems*2

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,connectg(1:totelems*2),dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

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

    ! Catchment ID
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/chi"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,chiID,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Node global ID
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/globID"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,drainers,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Strahler stream order
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/strahler"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,stNb,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(nodes,connect,facc,stNb,drainers,cID,connectg,chiID)

    return

  end subroutine drainage_hdf5
  ! =====================================================================================
  subroutine drainage_xmf(iter)

    type(xmlf_t)::xf

    integer::iter,totnodes,totelems,k,ierr

    character(len=128)::str,file,filename,filename1,filename2
    character(len=128)::stg,filename3,filename4,filename5

    call drainage_hdf5(iter,totelems)

    if(.not.allocated(doutnode)) allocate(doutnode(npets),doutelem(npets))
    call mpi_gather(drainOde,1,mpi_integer,doutnode,1,mpi_integer,0,badlands_world,ierr)
    call mpi_gather(totelems,1,mpi_integer,doutelem,1,mpi_integer,0,badlands_world,ierr)

    if(pet_id==0)then
        fdspm='DrainageSPM'
        totnodes=localNodes

        file=''
        file=fdspm
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
            totnodes=doutnode(k)
            totelems=doutelem(k)
            filename=''
            filename=fdspm
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
            str=':/connectivity'
            call append_str(filename,str)
            str=':/vertices'
            call append_str(filename1,str)
            str=':/facc'
            call append_str(filename2,str)
            str=':/cID'
            call append_str(filename3,str)
            str=':/strahler'
            call append_str(filename4,str)
            str=':/chi'
            call append_str(filename5,str)

            ! Block begin
            call xml_NewElement(xf,"Grid")
            str='DrainBlock.'
            call append_zero(str,iter)
            stg='.p'
            call append_str(str,stg)
            call append_zero(str,k-1)
            call xml_AddAttribute(xf,"Name",trim(str))
            call xml_NewElement(xf,"Topology")
            call xml_AddAttribute(xf,"Type","Polyline")
            call xml_AddAttribute(xf,"NodesPerElement",2)
            call xml_AddAttribute(xf,"NumberOfElements",totelems)
            call xml_AddAttribute(xf,"BaseOffset","1")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"DataType","Int")
            str=' '
            call append_nb2(str,totelems)
            call append_nb2(str,2)
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

            ! Strahler stream order
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Strahler stream order")
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

            ! Catchment ID
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","cID")
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

            ! Chi parameter
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","chi")
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
            call xml_EndElement(xf,"Grid")
        enddo

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    return

  end subroutine drainage_xmf
  ! =====================================================================================
  subroutine visualise_drainage_changes(iter)

    ! Parameters Declaration
    type(xmlf_t)::xf

    integer::i,iter,it0
    character(len=128)::filename,str,fname

    call drainage_xmf(iter)

    if(pet_id==0)then
        filename='SPMdrainage_series.xdmf'
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
            fname=fdspm
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

   end subroutine visualise_drainage_changes
   ! =====================================================================================
end module outspm_drainage
