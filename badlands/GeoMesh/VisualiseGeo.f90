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
!       Filename:  SurfaceOut.f90
!
!    Description:  Implements the XdmF and HdF5 delaunay and voronoi output grids
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module outsurf

  use hdf5
  use parameters
  use topology
  use FoX_wxml

  implicit none

  character(len=128)::fdelaunay,fvoronoi

contains

  ! =====================================================================================
  subroutine delaunay_hdf5

    logical::compression

    integer::id,i,rank,iter,totnodes,totelems
    integer,dimension(:),allocatable::connect

    real(kind=8),dimension(:),allocatable::nodes

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    iter=0

    fdelaunay='delaunay.h5'

    file=''
    file=fdelaunay
    call noblnk(file)
    call addpath1(file)
    totnodes=dnodes
    totelems=delemo

    allocate(nodes(3*totnodes))
    allocate(connect(3*totelems))

    ! Create nodes arrays
    id=1
    do i=1,totnodes
       nodes(id)=tcoordX(i)
       nodes(id+1)=tcoordY(i)
       nodes(id+2)=tcoordZ(i)
       id=id+3
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f( h5p_file_access_f, plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    ! Surface - connectivity
    id=1
    do i=1,delem
       if(elemtmask(i)==0)then
          connect(id)=delmt(i,1)
          connect(id+1)=delmt(i,2)
          connect(id+2)=delmt(i,3)
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

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(nodes,connect)

    return

  end subroutine delaunay_hdf5
  ! =====================================================================================
  subroutine delaunay_xmf

    type(xmlf_t)::xf

    integer::iter,totnodes,totelems

    character(len=128)::str,filename,filename1,filename2,file

    call delaunay_hdf5

    iter=0

    fdelaunay='delaunay'
    totnodes=dnodes
    totelems=delemo

    file=''
    file=fdelaunay
    call noblnk(file)
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
    call xml_AddAttribute(xf,"Value",0.0)
    call xml_EndElement(xf,"Time")

    filename=''
    filename=fdelaunay
    call noblnk(filename)
    str='.h5'
    call append_str(filename,str)
    filename1=filename
    filename2=filename
    str=':/connectivity'
    call append_str(filename,str)
    str=':/vertices'
    call append_str(filename1,str)

    ! Block begin
    call xml_NewElement(xf,"Grid")
    str='Block_delaunay.'
    call append_zero(str,iter)
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
    call xml_EndElement(xf,"Grid")

    ! Footer
    call xml_EndElement(xf,"Grid")
    call xml_EndElement(xf,"Domain")
    call xml_EndElement(xf,"Xdmf")
    call xml_Close(xf)

    return

  end subroutine delaunay_xmf
  ! =====================================================================================
  subroutine voronoi_hdf5

    logical :: compression

    ! Parameters Declaration
    integer::id,i,rank,iter,totnodes

    real(kind=8),dimension(:),allocatable::nodes

    character(len=128)::text,stg,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    iter=0

    fvoronoi='voronoi'

    file=''
    file=fvoronoi
    call noblnk(file)

    stg='.h5'
    call append_str(file,stg)
    call addpath1(file)
    totnodes=vnodes

    allocate(nodes(totnodes*3))

    ! Create nodes arrays
    id=1
    do i=1,vnodes
       nodes(id)=real(vcoordX(i))
       nodes(id+1)=real(vcoordY(i))
       nodes(id+2)=0.0
       id=id+3
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f(h5p_file_access_f,plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    ! The Coordinates - vertices
    dims(1)=3
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/vertices"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc)
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

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(nodes)

    return

  end subroutine voronoi_hdf5
  ! =====================================================================================
  subroutine voronoi_xmf

    type(xmlf_t)::xf

    integer::iter,totnodes,totelems,i,k

    character(len=128)::str,filename,filename1,file

    call voronoi_hdf5

    iter=0

    fvoronoi='voronoi'
    totnodes=vnodes
    totelems=vcellIN

    file=''
    file=fvoronoi
    call noblnk(file)
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
    call xml_AddAttribute(xf,"Value",0.0 )
    call xml_EndElement(xf,"Time")

    filename=''
    filename=fvoronoi
    call noblnk(filename)
    str = '.h5'
    call append_str(filename,str)
    filename1=filename
    str=':/vertices'
    call append_str(filename1,str)

    ! Block begin
    call xml_NewElement(xf,"Grid")
    str='Block_voronoi.'
    call append_zero(str,iter)
    call xml_AddAttribute(xf,"Name",trim(str))
    call xml_NewElement(xf,"Topology")
    call xml_AddAttribute(xf,"Type","Mixed")
    call xml_AddAttribute(xf,"NumberOfElements",totelems)
    call xml_NewElement(xf,"DataItem")
    call xml_AddAttribute(xf,"Format","XML")
    call xml_AddAttribute(xf,"DataType","Int")
    call xml_AddAttribute(xf,"Dimensions",velemIN+2*vcellIN)
    call xml_AddNewLine(xf)

    do i=1,dnodes
       if(voronoiCell(i)%border==0.and.voronoiCell(i)%vertexNb>0)then
          str='3 '
          call append_nb2(str,voronoiCell(i)%vertexNb)
          do k=1,voronoiCell(i)%vertexNb
             call append_nb2(str,voronoiCell(i)%vertexID(k)-1)
          enddo
          call xml_AddCharacters(xf,trim(str))
          call xml_AddNewLine(xf)
       endif
    enddo
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

    call xml_EndElement(xf,"Grid")

    ! Footer
    call xml_EndElement(xf,"Grid")
    call xml_EndElement(xf,"Domain")
    call xml_EndElement(xf,"Xdmf")
    call xml_Close(xf)

    return

  end subroutine voronoi_xmf
  ! =====================================================================================
end module outsurf
