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
!       Filename:  MeshPartition.f90
!
!    Description:  Implements the partitioning of unstructured delaunay grid & regular grid
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module meshPartition

  use parallel
  use parameters
  use topology
  use zoltan

  implicit none

  ! Data Recursive Coordinate Bisection (RCB)
  integer::numGlobObjs,numLocObjs
  integer,dimension(:),allocatable::zolt_part
  integer(ZOLTAN_INT),dimension(:),allocatable::zGID
  real,dimension(:),allocatable::zolX,zolY

  ! Zoltan data to store
  logical::changes
  integer(Zoltan_INT)::numGidEntries,numLidEntries
  integer(Zoltan_INT)::numImport,numExport
  integer(Zoltan_INT),pointer,dimension(:)::importGlobalGids,exportGlobalGids
  integer(Zoltan_INT),pointer,dimension(:)::importLocalGids,exportLocalGids
  integer(Zoltan_INT),pointer,dimension(:)::importProcs,exportProcs
  integer(Zoltan_INT),pointer,dimension(:)::importToPart,exportToPart

contains

  ! =====================================================================================
  subroutine StructureGridPart

    logical::column_partition

    integer::k,e,p,id,extra,partNb,partelems

    integer,dimension(npets)::partDim
    integer,dimension(nbelm)::elementPet
    integer,dimension(bnbnodes)::nodeDefine

    real(kind=8)::elem_center
    real(kind=8),dimension(npets)::maxCoord

    column_partition=.false.
    if(nx>ny)column_partition=.true.

    if(column_partition)then
      partelems=int((nx+1)/npets)
      extra=mod(nx+1,npets)
    else
      partelems=int((ny+1)/npets)
      extra=mod(ny+1,npets)
    endif

    partDim=0
    elementPet=-1
    if(pet_id==0)then
      ! Get the number of rows/columns for each PET
      do k=0,npets-1
        partNb=partelems
        if(k<extra)partNb=partNb+1
        partDim(k+1)=partNb+1
        if(k==0)then
          if(.not.column_partition)then
            maxCoord(k+1)=miny-dx+partNb*dx
          else
            maxCoord(k+1)=minx-dx+partNb*dx
          endif
        else
          maxCoord(k+1)=maxCoord(k)+partNb*dx
        endif
      enddo

      ! Get each element processor ID
      do e=1,nbelm
        loopface:do k=1,npets
          elem_center=rcoordX(relmt(e,1))+dx/2.0
          if(.not.column_partition) elem_center=rcoordY(relmt(e,1))+dx/2.0
          if(elem_center<maxCoord(k))then
            elementPet(e)=k-1
            exit loopface
          endif
        enddo loopface
      enddo
    endif
    call mpi_bcast(partDim,npets,mpi_integer,0,badlands_world,rc)
    call mpi_bcast(maxCoord,npets,mpi_double_precision,0,badlands_world,rc)
    call mpi_bcast(elementPet,nbelm,mpi_integer,0,badlands_world,rc)

    if(.not.column_partition)then
      spartN=partDim(pet_id+1)*(nx+2)
    else
      spartN=partDim(pet_id+1)*(ny+2)
    endif

    spartE=0
    do e=1,nbelm
      if(elementPet(e)==pet_id)spartE=spartE+1
    enddo

    if(allocated(snodeID)) deallocate(snodeID)
    if(allocated(snodeLID)) deallocate(snodeLID)
    allocate(snodeID(spartN),snodeLID(bnbnodes))
    id=0
    snodeLID=-1
    snodeID=-1
    do p=1,bnbnodes
      elem_center=rcoordX(p)
      if(.not.column_partition)elem_center=rcoordY(p)
      if(pet_id>0)then
        if(elem_center>=maxCoord(pet_id).and.elem_center<=maxCoord(pet_id+1))then
          id=id+1
          snodeID(id)=p
          snodeLID(p)=id
        endif
      endif
      if(pet_id==0.and.elem_center<=maxCoord(pet_id+1))then
        id=id+1
        snodeID(id)=p
        snodeLID(p)=id
      endif
    enddo

    if(allocated(selemID)) deallocate(selemID)
    if(allocated(sownedID)) deallocate(sownedID)
    allocate(sownedID(spartN),selemID(spartE))
    selemID=-1
    sownedID=pet_id
    nodeDefine=-1
    k=0
    id=0
    do e=1,nbelm
      if(elementPet(e)==pet_id)then
        k=k+1
        selemID(k)=e
      endif
    enddo

    if(pet_id>0)then
      do p=1,spartN
        elem_center=rcoordX(snodeID(p))
        if(.not.column_partition)elem_center=rcoordY(snodeID(p))
        if(maxCoord(pet_id)==elem_center)sownedID(p)=pet_id-1
      enddo
    endif

    sOwnedNode=0
    do p=1,spartN
      if(sownedID(p)==pet_id) sOwnedNode=sOwnedNode+1
    enddo

    return

  end subroutine StructureGridPart
  ! =====================================================================================
  subroutine UnstructureGridPart

    integer(Zoltan_INT)::ierr
    real(Zoltan_FLOAT)::version

    ierr=Zoltan_Initialize(version)
    allocate(zolt_part(delem))

    call allocateObjects
    call RCBpartition()
    call define_partition()

    deallocate(zGID)
    deallocate(zolX,zolY)
    call zoltanCleanUp()
    deallocate(zolt_part)

    return

  end subroutine UnstructureGridPart
  ! =====================================================================================
  subroutine allocateObjects

    integer::i,currIndx

    numGlobObjs=delem
    numLocObjs=0

    ! Round robin initial distribution
    do i=1,numGlobObjs
       if(mod(i-1,npets)==pet_id) numLocObjs=numLocObjs+1
    enddo

    ! Allocate data for round robin distribution grid
    allocate(zGID(numLocObjs))
    allocate(zolX(numLocObjs))
    allocate(zolY(numLocObjs))

    ! Fill data for round robin distribution grid
    currIndx=1
    do i=1,numGlobObjs
       ! assumes gids start at 1, gives round robin initial distribution
       if(mod(i-1,npets)==pet_id)then
          zGID(currIndx)=i
          zolX(currIndx)=real(dcentroid(i,1))
          zolY(currIndx)=real(dcentroid(i,2))
          currIndx=currIndx+1
       endif
    enddo

  end subroutine allocateObjects
  ! =====================================================================================

  subroutine define_partition

    integer::partAssign(numGlobObjs)
    integer::allPartAssign(numGlobObjs)
    integer::parts(numLocObjs)
    integer::i,p,id
    integer,dimension(:),allocatable::upartIDpts,upartIDelem,uIDpts

    do i=1,numLocObjs,1
       parts(i)=pet_id
    enddo

    do i=1,numExport
       parts(exportLocalGids(i))=exportToPart(i)
    enddo

    partAssign=0
    allPartAssign=0
    zolt_part=0

    do i=1,numLocObjs
       partAssign(zGID(i))=parts(i)
    enddo
    call mpi_allreduce(partAssign,allPartAssign,numGlobObjs,mpi_integer,mpi_max,badlands_world,rc)

    if(pet_id==0)then
       do i=1,numGlobObjs
          zolt_part(i)=allPartAssign(i)
       enddo
    endif
    call mpi_bcast(zolt_part,numGlobObjs,mpi_integer,0,badlands_world,rc)

    ! Defined the unstructured mesh partition parameters
    upartE=0
    upartN=0
    if(allocated(upartIDelem)) deallocate(upartIDelem)
    if(allocated(upartIDpts)) deallocate(upartIDpts)
    if(allocated(uIDpts)) deallocate(uIDpts)
    if(allocated(uownEID)) deallocate(uownEID)
    allocate(upartIDpts(dnodes),uIDpts(dnodes),upartIDelem(delem),uownEID(delem))
    upartIDpts=-1
    uIDpts=-1
    do i=1,delem
      upartIDelem(i)=zolt_part(i)
      if(zolt_part(i)==pet_id) upartE=upartE+1
      do p=1,3
        if(zolt_part(i)==pet_id)then
          if(uIDpts(delmt(i,p))==-1)then
            uIDpts(delmt(i,p))=pet_id
            upartN=upartN+1
          endif
        endif
        if(upartIDpts(delmt(i,p))<0)then
          upartIDpts(delmt(i,p))=zolt_part(i)
        endif
      enddo
    enddo

    ! Get the partition elements
    if(allocated(uelemID)) deallocate(uelemID)
    allocate(uelemID(upartE))
    p=0
    uownEID=-1
    do i =1,delem
      if(upartIDelem(i)==pet_id)then
        p=p+1
        uelemID(p)=i
        uownEID(i)=pet_id
      endif
    enddo

    ! Get the partition nodes and owned processor ID
    if(allocated(unodeID)) deallocate(unodeID)
    if(allocated(unodeID)) deallocate(unodeID)
    if(allocated(unodeLID)) deallocate(unodeLID)
    if(allocated(uownedID)) deallocate(uownedID)
    allocate(unodeID(upartN),uownedID(upartN),unodeLID(dnodes))
    unodeLID=-1
    p=0
    id=0
    do i=1,dnodes
      if(uIDpts(i)==pet_id)then
        p=p+1
        unodeID(p)=i
        unodeLID(i)=p
        uownedID(p)=upartIDpts(i)
        if(uownedID(p)==pet_id)id=id+1
      endif
    enddo
    uOwnedNode=id

    deallocate(uIDpts,upartIDpts,upartIDelem)

  end subroutine define_partition
  ! =====================================================================================

  subroutine RCBPartition

    ! Local variables
    type(Zoltan_Struct),pointer::zz_obj
    integer(Zoltan_INT)::ierr

    ! Body of subroutine
    nullify(zz_obj)

    zz_obj=>Zoltan_Create(badlands_world)

    ! General Zoltan Parameters
    ierr=Zoltan_Set_Param(zz_obj,"DEBUG_LEVEL","0")
    ierr=Zoltan_Set_Param(zz_obj,"LB_METHOD","RCB")

    ! Register query functions
    ierr=Zoltan_Set_Fn(zz_obj,ZOLTAN_NUM_OBJ_FN_TYPE,zoltNumObjs)
    ierr=Zoltan_Set_Fn(zz_obj,ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
    ierr=Zoltan_Set_Fn(zz_obj,ZOLTAN_NUM_GEOM_FN_TYPE,zoltNumGeom)
    ierr=Zoltan_Set_Fn(zz_obj,ZOLTAN_GEOM_FN_TYPE,zoltGeom)

    ierr=Zoltan_LB_Partition(zz_obj,changes,numGidEntries,numLidEntries, &
         numImport,importGlobalGids,importLocalGids,importProcs,importToPart, &
         numExport,exportGlobalGids,exportLocalGids,exportProcs,exportToPart)

    call Zoltan_Destroy(zz_obj)

  end subroutine RCBPartition
  ! =====================================================================================
  subroutine zoltanCleanup()

    integer::ierr

    ierr=Zoltan_LB_Free_Part(importGlobalGids,importLocalGids,importProcs,importToPart)
    ierr=Zoltan_LB_Free_Part(exportGlobalGids,exportLocalGids,exportProcs,exportToPart)

  end subroutine zoltanCleanup
  ! =====================================================================================
  integer function zoltNumObjs(data,ierr)

    ! Local declarations
    integer(Zoltan_INT),intent(in)::data(*)
    integer(Zoltan_INT),intent(out)::ierr

    zoltNumObjs=numLocObjs
    ierr=ZOLTAN_OK

  end function zoltNumObjs
  ! =====================================================================================
  subroutine zoltGetObjs(data,num_gid_entries,num_lid_entries,global_ids, &
       local_ids,wgt_dim,obj_wgts,ierr)

    integer(Zoltan_INT),intent(in)::data(*)
    integer(Zoltan_INT),intent(in)::num_gid_entries
    integer(Zoltan_INT),intent(in)::num_lid_entries
    integer(Zoltan_INT),intent(out)::global_ids(*)
    integer(Zoltan_INT),intent(out)::local_ids(*)
    integer(Zoltan_INT),intent(in)::wgt_dim
    real(Zoltan_FLOAT),intent(out)::obj_wgts(*)
    integer(Zoltan_INT),intent(out)::ierr

    integer::i

    do i=1,numLocObjs
       global_ids(i)=zGID(i)
       local_ids(i)=i
    enddo

    ierr=ZOLTAN_OK

  end subroutine zoltGetObjs
  ! =====================================================================================
  integer function zoltNumGeom(data,ierr)

    integer(Zoltan_INT),intent(in)::data(*)
    integer(Zoltan_INT)::ierr

    zoltNumGeom=2
    ierr=ZOLTAN_OK

  end function zoltNumGeom
  ! =====================================================================================

  subroutine zoltGeom(data,num_gid_entries,num_lid_entries,global_id, &
       local_id,geom_vec,ierr)

    integer(Zoltan_INT),intent(in)::data(*)
    integer(Zoltan_INT),intent(in)::num_gid_entries
    integer(Zoltan_INT),intent(in)::num_lid_entries
    integer(Zoltan_INT),intent(in)::global_id
    integer(Zoltan_INT),intent(in)::local_id
    real(Zoltan_DOUBLE),intent(out)::geom_vec(*)
    integer(Zoltan_INT),intent(out)::ierr

    geom_vec(1)=zolX(local_id)
    geom_vec(2)=zolY(local_id)

    ierr=ZOLTAN_OK

  end subroutine zoltGeom
  ! =====================================================================================

end module meshPartition
! =====================================================================================
