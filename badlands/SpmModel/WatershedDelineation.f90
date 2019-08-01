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
!       Filename:  WatershedDelineation.f90
!
!    Description:  Sub-basin partitioning based upon channel network junctions
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module watershed

  use parallel
  use topology
  use parameters
  use hydroUtil
  use hydrology

  implicit none

  integer,dimension(:),allocatable::rcvNb
  integer,dimension(:,:),allocatable::rcvIDs

  integer::junctionsNb,jctNb
  integer,dimension(:),allocatable::junctionIDs,ljunctionIDs
  integer,dimension(:),allocatable::partition

  integer,dimension(:),allocatable::options
  integer,pointer::vsize,adjwgt
  real(kind=8),pointer::tpwgts,ubvec

contains

  ! =====================================================================================
  subroutine compute_subcatchment

    integer::id,k,n,p,s,maxs,jcts(npets),disp(npets),disps(npets+1)
    ! real(kind=8)::act,aire
    real(kind=8)::pp,dist
    if(.not.allocated(strahler)) allocate(strahler(dnodes))
    if(.not.allocated(subcatchmentID)) allocate(subcatchmentID(dnodes))
    if(.not.allocated(junctionIDs)) allocate(junctionIDs(dnodes))
    if(.not.allocated(ljunctionIDs)) allocate(ljunctionIDs(dnodes))

    if(.not.allocated(rcvNb)) allocate(rcvNb(dnodes))
    if(.not.allocated(rcvIDs)) allocate(rcvIDs(dnodes,maxrcvs*2))
    rcvNb=0
    rcvIDs=0

    ! aire=del_area
    ! if(refineNb>0)then
    !   aire=refine_grid(1)%area
    ! endif

    ! Compute cumulative discharge
    do id=partStack(pet_id+1),1,-1
      k=lstackOrder(id)
      p=receivers(k)
      ! Receivers
      rcvNb(p)=rcvNb(p)+1
      rcvIDs(p,rcvNb(p))=k
      if(p/=k) discharge(p)=discharge(p)+discharge(k)
    enddo

    call mpi_allreduce(mpi_in_place,discharge,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

    ! Create global catchment IDs and Chi coefficient
    s=0
    pp=spl_m/spl_n
    do id=1,dnodes
      k=stackOrder(id)
      p=receivers(k)
      if(p==k) s=s+1
      bsID(k)=s
      if(discharge(p)>0..and.discharge(k)>0.and.p/=k)then
        dist=sqrt((tcoordX(k)-tcoordX(p))**2.+(tcoordY(k)-tcoordY(p))**2.)
        chi(k)=chi(p)+0.5*((1./(discharge(p)))**pp+(1./(discharge(k)))**pp)*dist
      endif
    enddo


    ! Compute Strahler stream order
    strahler=0
    jctNb=0
    ljunctionIDs=0
    ! act=aire*accu_thres
    do id=partStack(pet_id+1),1,-1
      maxs=0
      k=lstackOrder(id)
      ! outlet
      n=0
      do p=1,rcvNb(k)
        s=strahler(rcvIDs(k,p))
        if(s>=maxs)then
          if(s==maxs)then
            n=n+1
          else
            n=1
          endif
          maxs=s
        endif
      enddo
      if(n==1)strahler(k)=maxs
      if(n>1)strahler(k)=maxs+1

      if(receivers(k)==k)then
        jctNb=jctNb+1
        ljunctionIDs(jctNb)=k
        ! act=aire*accu_thres !accu_thres
      else
        if(n>1.and.strahler(k)>1.and.mod(strahler(k),2)==0)then
        ! if(discharge(k)>act)then
          ! act=act+aire*accu_thres
          jctNb=jctNb+1
          ljunctionIDs(jctNb)=k
        endif
      endif
    enddo

    call mpi_allreduce(mpi_in_place,strahler,dnodes,mpi_integer,mpi_max,badlands_world,rc)
    call mpi_allreduce(jctNb,junctionsNb,1,mpi_integer,mpi_sum,badlands_world,rc)
    call mpi_allgather(jctNb,1,mpi_integer,jcts,1,mpi_integer,badlands_world,rc)
    disp=0
    disps=0
    do p=1,npets
      if(p<npets) disp(p+1)=disp(p)+jcts(p)
      disps(p+1)=disps(p)+jcts(p)
    enddo
    junctionIDs=0
    call mpi_allgatherv(ljunctionIDs,jctNb,mpi_integer,junctionIDs,jcts,disp,mpi_integer,badlands_world,rc)

    ! Count nodes per subcatchment
    if(allocated(subcatchNb)) deallocate(subcatchNb)
    allocate(subcatchNb(junctionsNb))
    subcatchNb=0
    subcatchmentID=-1
    id=0
    do n=disps(pet_id+1)+1,disps(pet_id+2)
      id=id+1
      s=ljunctionIDs(id)
      subcatchmentID(s)=n
      subcatchNb(n)=subcatchNb(n)+1
      p=addtosubcatch(s,n)
    enddo

    call mpi_allreduce(mpi_in_place,junctionIDs,dnodes,mpi_integer,mpi_max,badlands_world,rc)
    call mpi_allreduce(mpi_in_place,subcatchNb,junctionsNb,mpi_integer,mpi_max,badlands_world,rc)
    call mpi_allreduce(mpi_in_place,subcatchmentID,dnodes,mpi_integer,mpi_max,badlands_world,rc)

    ! Perform subcatchment load-balancing and partitioning
    if(pet_id==0)call metis_loadbalancing

    return

  end subroutine compute_subcatchment
  ! =====================================================================================
  recursive function addtosubcatch(k,catchID) result(success)

    integer::n,success,k,catchID,s

    success=1

    do n=1,rcvNb(k)
      s=rcvIDs(k,n)
      if(subcatchmentID(s)==-1)then
        subcatchmentID(s)=catchID
        subcatchNb(catchID)=subcatchNb(catchID)+1
        success=addtosubcatch(s,catchID)
      endif
    enddo

    success=0

  end function addtosubcatch
  ! =====================================================================================
  subroutine metis_loadbalancing

    integer::l,s,id,k,p,objval

    integer,dimension(:),allocatable::vwgt,CRSadj,CRSadjncy,connectsNb
    integer,dimension(:,:),allocatable::connects

    if(allocated(partition)) deallocate(partition)
    allocate(partition(junctionsNb+1))

    if(allocated(vwgt)) deallocate(vwgt,CRSadj,CRSadjncy)
    if(allocated(connects)) deallocate(connects,connectsNb)
    allocate(vwgt(junctionsNb+1))
    allocate(CRSadj(junctionsNb+2),CRSadjncy(2*(junctionsNb)))
    allocate(connectsNb(junctionsNb+1),connects(junctionsNb+1,junctionsNb+1))
    connectsNb=0

    do s=1,junctionsNb+1
      if(s<=junctionsNb)then
        id=junctionIDs(s)
        k=receivers(id)
        l=subcatchmentID(k)
        vwgt(s)=subcatchNb(s)
      else
        l=s
        vwgt(s)=1
      endif
      ! Declare tree connection parameters
      if(l/=s)then
        connectsNb(s)=connectsNb(s)+1
        connectsNb(l)=connectsNb(l)+1
        connects(s,connectsNb(s))=l
        connects(l,connectsNb(l))=s
      ! Link to the top root
      elseif(s<=junctionsNb)then
        l=junctionsNb+1
        connectsNb(s)=connectsNb(s)+1
        connectsNb(l)=connectsNb(l)+1
        connects(s,connectsNb(s))=l
        connects(l,connectsNb(l))=s
      endif
    enddo

    ! Define the metis graph data structure CRS
    CRSadj=1
    do s=2,junctionsNb+2
      CRSadj(s)=CRSadj(s-1)+connectsNb(s-1)
    enddo

    CRSadjncy=0
    k=1
    do p=1,junctionsNb+1
      do l=1,connectsNb(p)
        CRSadjncy(k)=connects(p,l)
        k=k+1
      enddo
    enddo

    ! Metis option numbering
    allocate(options(0:40))
    call METIS_SetDefaultOptions(options)
    options(17)=1

    ! Nullify pointers
    vsize=>null()
    adjwgt=>null()
    tpwgts=>null()
    ubvec=>null()

    ! K-way metis partitioning
    partition=1
    if(npets>1) call METIS_PartGraphKway(junctionsNb+1,1,CRSadj,CRSadjncy,vwgt,vsize,adjwgt,npets,tpwgts,ubvec,options,objval,partition)
    deallocate(options)

    ! Cleaning process
    nullify(vsize,adjwgt,tpwgts,ubvec)

  end subroutine metis_loadbalancing
  ! =====================================================================================
  subroutine bcast_loadbalancing

    integer::s,id,k,p

    ! Broadcast partitioning
    if(pet_id/=0)then
      if(allocated(partition)) deallocate(partition)
      allocate(partition(junctionsNb+1))
    endif
    call mpi_bcast(partition,junctionsNb+1,mpi_integer,0,badlands_world,rc)

    if(.not.allocated(subcatchmentProc)) allocate(subcatchmentProc(dnodes))
    subcatchmentProc=-1

    ! Define local and inter-processor communications
    if(.not.allocated(sendprocID)) allocate(sendprocID(dnodes))
    sendprocID=-1
    if(.not.allocated(rcvprocNb)) allocate(rcvprocNb(dnodes))
    rcvprocNb=0
    if(.not.allocated(rcvprocID)) allocate(rcvprocID(dnodes,maxrcvs))
    rcvprocID=-1
    if(.not.allocated(rcvsendID)) allocate(rcvsendID(dnodes,maxrcvs))
    rcvsendID=-1
    if(.not.allocated(localNodesGID)) allocate(localNodesGID(dnodes))

    ! Define for each partition local nodes and their global IDs
    s=1
    drainOde=0
    localNodes=0
    do k=1,dnodes
      p=stackOrder(k)
      id=receivers(p)
      subcatchmentProc(p)=partition(subcatchmentID(p))-1
      subcatchmentProc(id)=partition(subcatchmentID(id))-1
      subcatchmentProc(k)=partition(subcatchmentID(k))-1
      if(pet_id==subcatchmentProc(k))localNodes=localNodes+1
      ! Inter-partition nodes variables
      if(subcatchmentProc(p)/=subcatchmentProc(id))then
        sendprocID(p)=subcatchmentProc(id)
        rcvprocNb(id)=rcvprocNb(id)+1
        rcvprocID(id,rcvprocNb(id))=subcatchmentProc(p)
        rcvsendID(id,rcvprocNb(id))=p
        drainOde=drainOde+1
      endif
      if(pet_id==subcatchmentProc(p))then
        localNodesGID(s)=k
        s=s+1
      endif
    enddo

    drainOde=drainOde+localNodes


  end subroutine bcast_loadbalancing
  ! =====================================================================================
end module watershed
