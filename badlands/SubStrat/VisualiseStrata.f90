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
!       Filename:  VisualiseStrata.f90
!
!    Description:  Visualise stratigraphic grid
!
!        Version:  1.0
!        Created:  03/06/15 09:05:27
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module out_stratal

  use hdf5
  use parallel
  use parameters
  use FoX_wxml
  use stratal_class

  implicit none

  character(len=128)::fstrata

contains

  ! =====================================================================================
  subroutine strata_hdf5(iter)


    logical::compression

    integer::id,idt,gid,k,i,ks,p,iter,totnodes,ierr,rank,totelems

    integer,dimension(:),allocatable::layID,connect

    real(kind=8),dimension(upartN)::lzh
    real(kind=8),dimension(:),allocatable::nodes,thick,mz
    real(kind=8),dimension(:,:),allocatable::prop

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    fstrata='StratalField'

    file=''
    file=fstrata
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

    totnodes=upartN*layNb
    totelems=delemoo*(layNb-1)
    allocate(nodes(3*totnodes))
    allocate(prop(totnodes,totgrn))
    allocate(thick(totnodes))
    allocate(mz(totnodes))
    allocate(layID(totnodes))

    allocate(connect(6*totelems))

    ! Create nodes arrays
    id=1
    idt=1
    do k=1,upartN
       gid=unodeID(k)
       nodes(id)=tcoordX(gid)
       nodes(id+1)=tcoordY(gid)
       nodes(id+2)=lay_base(k)
       lzh(k)=lay_base(k)
       thick(idt)=0.0_8
       prop(idt,1:totgrn)=0.0_8
       layID(idt)=0
       mz(idt)=0.0_8
       id=id+3
       idt=idt+1
    enddo

    do k=1,layNb-1
      do i=1,upartN
        gid=unodeID(i)
        nodes(id)=tcoordX(gid)
        nodes(id+1)=tcoordY(gid)
        thick(idt)=lay_thick(i,k)
        nodes(id+2)=lzh(i)+thick(idt)
        lzh(i)=nodes(id+2)
        layID(idt)=k
        mz(idt)=0.0_8
        if(lay_thick(i,k)>0.)then
            do ks=1,totgrn
                prop(idt,ks)=lay_sed(i,k,ks)/lay_thick(i,k)
                mz(idt)=mz(idt)+prop(idt,ks)*sediments(ks)%dia
            enddo
        else
            prop(idt,1:totgrn)=0.0
        endif
        id=id+3
        idt=idt+1
      enddo
    enddo

    ! Surface - connectivity
    id=1
    do k=1,layNb-1
      do p=1,delemoo
        i=delemID(p)
        if(elemtmask(i)==0)then
          connect(id)=unodeLID(delmt(i,1))+(k-1)*upartN
          connect(id+1)=unodeLID(delmt(i,2))+(k-1)*upartN
          connect(id+2)=unodeLID(delmt(i,3))+(k-1)*upartN
          connect(id+3)=connect(id)+upartN
          connect(id+4)=connect(id+1)+upartN
          connect(id+5)=connect(id+2)+upartN
          id=id+6
        endif
      enddo
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f(h5p_file_access_f,plist_id,rc)
    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    dims(1)=6
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
    dims(2)=totelems*6

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,connect,dims,rc)
    call h5pclose_f(plist_id,rc)
    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! The Coordinates - vertices
    dims(1)=3
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/vertices"
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)
    dims(1)=1
    dims(2)=totnodes*3

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,nodes,dims,ierr)
    call h5pclose_f(plist_id,ierr)
    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Thickness attribute
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/thick"
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,thick,dims,ierr)
    call h5pclose_f(plist_id,ierr)
    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Layer number type attribute
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/layerID"
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,ierr,plist_id)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,layID,dims,ierr)
    call h5pclose_f(plist_id,ierr)
    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Mean grainsize attribute
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/mgz"
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)
    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,mz,dims,ierr)
    call h5pclose_f(plist_id,ierr)
    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Sediment proportion attribute
    do ks=1,totgrn
      dims(1)=1
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/prop"
      call append_nb(text,ks)
      ! Create property list for collective dataset write
      call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
      call h5pset_deflate_f(plist_id,9,ierr)
      call h5pset_chunk_f(plist_id,rank,dims,ierr)
      ! Create the dataset with default properties
      call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)
      ! Write the dataset collectively
      call h5dwrite_f(dset_id,h5t_native_double,prop(1:totnodes,ks),dims,ierr)
      call h5pclose_f(plist_id,ierr)
      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)
    enddo

    ! Close the file.
    call h5fclose_f(file_id,rc)
    ! Close interface
    call h5close_f(rc)

    deallocate(nodes,prop,thick,mz,layID,connect)

    return

  end subroutine strata_hdf5
  ! =====================================================================================
  subroutine strata_xmf(iter)

    type(xmlf_t)::xf
    integer::iter,totnodes,totelems,k,ks

    character(len=128)::str,stg,filename,filename1,filename2,file,filename3
    character(len=128)::filename4,filename5,filename6,txt

    call strata_hdf5(iter)

    if(pet_id==0)then
        fstrata='StratalField'
        file=''
        file=fstrata
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
            totnodes=outnode(k)*layNb
            totelems=outelem(k)*(layNb-1)
            filename=''
            filename=fstrata
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
            str=':/connectivity'
            call append_str(filename,str)
            str=':/vertices'
            call append_str(filename1,str)
            str=':/thick'
            call append_str(filename2,str)
            str=':/mgz'
            call append_str(filename3,str)
            str=':/layerID'
            call append_str(filename4,str)
            str=':/prop'
            call append_str(filename5,str)

            ! Block begin
!             call xml_NewElement(xf,"Grid")
!             str='StratBlock.'
!             call append_zero(str,iter)
!             stg='.p'
!             call append_str(str,stg)
!             call append_zero(str,k-1)
!             call xml_AddAttribute(xf,"Name",trim(str))
!             call xml_NewElement(xf,"Topology")
!             call xml_AddAttribute(xf,"TopologyType","3DSMesh")
!             str=' '
!             call append_nb2(str,nbZ)
!             call append_nb2(str,nbY(k))
!             call append_nb2(str,nbX(k))
!             call xml_AddAttribute(xf,"Dimensions",trim(str))
!             call xml_EndElement(xf,"Topology")

            ! Block begin
            call xml_NewElement(xf,"Grid")
            str='SurfBlock.'
            call append_zero(str,iter)
            stg='.p'
            call append_str(str,stg)
            call append_zero(str,k-1)
            call xml_AddAttribute(xf,"Name",trim(str))
            call xml_NewElement(xf,"Topology")
            call xml_AddAttribute(xf,"Type","Wedge")
            call xml_AddAttribute(xf,"NumberOfElements",totelems)
            call xml_AddAttribute(xf,"BaseOffset","1")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"DataType","Int")
            str=' '
            call append_nb2(str,totelems)
            call append_nb2(str,6)
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

            ! Layer thickness
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Layer Thickness")
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

            ! Mean grain size
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Mean Grainsize")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","4")
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename3))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Layer Nb
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Node")
            call xml_AddAttribute(xf,"Name","Layer Index")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Int")
            call xml_AddAttribute(xf,"Precision","4")
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename4))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Sediment proportion
            do ks=1,totgrn
                filename6=filename5
                call append_nb(filename6,ks)
                txt='mat'
                call append_nb(txt,ks)
                call xml_NewElement(xf,"Attribute")
                call xml_AddAttribute(xf,"Type","Scalar")
                call xml_AddAttribute(xf,"Center","Node")
                call xml_AddAttribute(xf,"Name",trim(txt))
                call xml_NewElement(xf,"DataItem")
                call xml_AddAttribute(xf,"Format","HDF")
                call xml_AddAttribute(xf,"NumberType","Float")
                call xml_AddAttribute(xf,"Precision","4")
                call xml_AddAttribute(xf,"Dimensions",trim(str))
                call xml_AddCharacters(xf,trim(filename6))
                call xml_EndElement(xf,"DataItem")
                call xml_EndElement(xf,"Attribute")
            enddo
            call xml_EndElement(xf,"Grid")
        enddo

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    return

  end subroutine strata_xmf
  ! =====================================================================================
  subroutine visualise_strata(iter)

    ! Parameters Declaration
    type(xmlf_t)::xf

    integer::i,iter,it0
    character(len=128)::filename,str,fname

    call strata_xmf(iter)

    if(pet_id==0)then
        filename='Stratal_series.xdmf'
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
            fname=fstrata
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

  end subroutine visualise_strata
  ! =====================================================================================
end module out_stratal
