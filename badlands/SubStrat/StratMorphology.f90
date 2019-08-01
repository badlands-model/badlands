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
!       Filename:  StratMorphology.f90
!
!    Description:  Defines morphological formalism used in BADLANDS for stratigraphic grid.
!
!        Version:  1.0
!        Created:  11/06/15 09:11:25
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module stratmorph

  use bilinear
  use parallel
  use topology
  use parameters
  use hydrology
  use hydroUtil
  use watershed
  use out_stratal
  use strata_evol
  use stratal_class
  use strat_hillslope
  use outspm_surface
  use outspm_drainage
  use external_forces

  implicit none

  integer::iter,idiff

  real(kind=8)::cpl_time,max_time
  real(kind=8),dimension(:),allocatable::nth
  real(kind=8),dimension(:,:),allocatable::nsed
!   real(kind=8)::time1,time2,tt1,tt2

contains

  ! =====================================================================================
  subroutine stratgeomorph

    integer::k,ks

    real(kind=8),dimension(dnodes)::tsed,gsed

    if(.not.allocated(nZ)) allocate(nZ(dnodes))
    if(.not.allocated(nth)) allocate(nth(dnodes))
    if(.not.allocated(newZ)) allocate(newZ(dnodes))
    if(.not.allocated(spmZ)) allocate(spmZ(dnodes))
    if(.not.allocated(spmH)) allocate(spmH(dnodes))
    if(.not.allocated(nsed)) allocate(nsed(dnodes,totgrn))
    if(.not.allocated(Qs_inS)) allocate(Qs_inS(dnodes,totgrn))
    if(.not.allocated(change_localS)) allocate(change_localS(dnodes,totgrn))

    ! Define paramaters
    if(simulation_time==time_start.or.update3d)then
      if(simulation_time==time_start) iter=0
      ! Build active layer
      call buildActiveLayer
      spmZ=tcoordZ
      CFL_diffusion=display_interval
      idiff=0
      if(simulation_time==time_start) time_step=0.
      if(simulation_time==time_start) time_display=time_start
      if(simulation_time==time_start) time_strata=time_start
      if(Cefficiency>0..or.stream_ero>0.) Tforce=1
    endif

    cpl_time=min(cpl1_time,cpl2_time)
    cpl_time=min(cpl_time,time_end)

    do while(simulation_time<cpl_time+0.001)

      if(.not.allocated(filldem)) allocate(filldem(dnodes))
      if(.not.allocated(watercell)) allocate(watercell(dnodes))

      ! Update borders
      call update_grid_borders

      ! Perform depressionless water filling algo
      call planchon_dem_fill_algorithm

      ! Find network tree based on Braun & Willet 2013
      call define_landscape_network

      ! Define subcathcment partitioning
      call compute_subcatchment

      ! Define load balancing
      call bcast_loadbalancing

      if(simulation_time==time_start.or.update3d) newZ=spmZ
      if(pet_id==0)print*,'Current time:',simulation_time

      ! Add stratigraphic layer
      if(simulation_time>=time_strata)then
        if(time_strata<=time_end) layNb=layNb+1
        time_strata=time_strata+time_layer
      endif

      ! Visualisation surface
      if(simulation_time>=time_display)then
        call visualise_surface_changes(iter)
        call visualise_drainage_changes(iter)
        call visualise_strata(iter)
        call mpi_barrier(badlands_world,rc)
        if(pet_id==0)print*,'Creating output: ',int(simulation_time)
        time_display=time_display+display_interval
        iter=iter+1
        newZ=spmZ
      endif

      ! Get time step size for hillslope process and stream power law
      call CFL_conditionS

      ! Geomorphological evolution
      call geomorphic_evolutionS

      ! Advance time
      simulation_time=simulation_time+time_step
      ! Update sea-level
      if(gsea%sealevel) call eustatism
      ! Apply displacement
      if(disp%event>0.and..not.disp3d)then
        call compute_vertical_displacement
        call compute_stratal_displacement
      endif

      ! Merge local geomorphic evolution
      call mpi_allreduce(nZ,spmZ,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

      ! Merge local stratigraphic evolution
      call mpi_allreduce(nth,alay_thick,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

      do ks=1,totgrn
        tsed(1:dnodes)=nsed(1:dnodes,ks)
        call mpi_allreduce(tsed,gsed,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
        alay_sed(1:dnodes,ks)=gsed(1:dnodes)
      enddo

      ! Update stratigraphic layer
      call update_stratigraphy_layer
      update3d=.false.
!       print*,time_step
!       stop

    enddo

    if(simulation_time>=time_end)then
      if(time_display<=time_end)then
        call visualise_surface_changes(iter)
        call visualise_drainage_changes(iter)
        call visualise_strata(iter)
      endif
      if(pet_id==0)print*,'simulation time: ',int(time_end)
    endif

    ! Update elevation for coupling component assuming infinite boundary
    do k=1,dnodes
      if(voronoiCell(k)%border==0) tcoordZ(k)=spmZ(k)
    enddo
    call update_grid_borders

  end subroutine stratgeomorph
  ! =====================================================================================
  subroutine CFL_conditionS

    integer::k,p,id,lid
    real(kind=8)::denom,distance,dh,ldt,dt

    ! Hillslope CFL conditions
    if(simulation_time==time_start.or.idiff==0.or.sediments(1)%diffa>0.)then
      call CFLdiffusionS
      CFL_diffusion=display_interval
      idiff=1
      if(CFL_diffusion==0.) time_step=display_interval
    else
      time_step=CFL_diffusion
    endif

    ! CFL factor for stream power law (bedrock incision)
    if(Cerodibility>0.)then
      do lid=1,localNodes
        p=localNodesGID(lid)
        k=stackOrder(p)
        if(voronoiCell(k)%border==0)then
          id=receivers(k)
          if(k/=id)then
            dh=spmZ(k)-spmZ(id)
            distance=sqrt((tcoordX(k)-tcoordX(id))**2.0+(tcoordY(k)-tcoordY(id))**2.0)
            if(spmZ(k)>spmZ(id).and.discharge(k)>0.)then
              denom=Cerodibility*discharge(k)**spl_m*(dh/distance)**(spl_n-1.0)
              time_step=min(time_step,distance/denom)
            endif
          endif
        endif
      enddo
    endif

    if(time_step>display_interval) time_step=display_interval
    if(time_step<force_time) time_step=force_time

    ! Get maximum time step for stability
    ldt=time_step
    call mpi_allreduce(ldt,dt,1,mpi_double_precision,mpi_min,badlands_world,rc)
    time_step=dt

    ! Check time-step in relation to coupling time
    if(simulation_time+time_step>cpl_time) time_step=cpl_time-simulation_time+1.e-4
    if(simulation_time+time_step>time_display) time_step=time_display-simulation_time+1.e-4
    if(simulation_time+time_step>time_strata) time_step=time_strata-simulation_time+1.e-4
    if(simulation_time+time_step>time_end) time_step=time_end-simulation_time+1.e-4

  end subroutine CFL_conditionS
  ! =====================================================================================
  subroutine geomorphic_evolutionS

    integer::k,id,rcv,lid,p,q,m,ks
    integer::stat(mpi_status_size),ierr,req(localNodes),r

    real(kind=8),dimension(totgrn)::LDL,SPL,Qs1,Qsr

    real(kind=8)::mtime,dt,maxtime,th
    real(kind=8)::distance,diffH,maxh

    max_time=time_step
    Qs_inS=0.0
    nZ=-1.e6
    nth=0.0
    nsed=0.0
    change_localS=-1.e6
    r=0

    do lid=localNodes,1,-1
      k=localNodesGID(lid)
      id=stackOrder(k)

      ! Receive child influx
      if(rcvprocNb(id)>0)then
        do p=1,rcvprocNb(id)
          call mpi_recv(Qsr,totgrn,mpi_double_precision,rcvprocID(id,p),rcvsendID(id,p),badlands_world,stat,ierr)
          do ks=1,totgrn
            Qs_inS(id,ks)=Qs_inS(id,ks)+Qsr(ks)
          enddo
        enddo
      endif

      if(voronoiCell(id)%btype<0)then

        ! Define local parameters
        rcv=receivers(id)
        distance=sqrt((tcoordX(id)-tcoordX(rcv))**2.0+(tcoordY(id)-tcoordY(rcv))**2.0)

        ! Creep processes
        call hillslope_fluxS(id,LDL,diffH)

        ! Stream Power Law (detachment-limited) - bedrock incision
        SPL=0.
        Qs1=0.
        call detachmentlimitedS(id,distance,diffH,SPL,Qs1)

        ! Detachment-limited condition
        do ks=1,totgrn
          Qs_inS(rcv,ks)=Qs_inS(rcv,ks)+Qs1(ks)
        enddo

        do ks=1,totgrn
          if(perosive==1.and.SPL(ks)>0.) SPL(ks)=0.
          change_localS(id,ks)=SPL(ks)+LDL(ks)
        enddo

      else
        change_localS(id,1:totgrn)=0.
      endif

      ! Send parent outflux
      if(sendprocID(id)>=0)then
        r=r+1
        call mpi_isend(Qs_inS(rcv,1:totgrn),totgrn,mpi_double_precision,sendprocID(id),id,badlands_world,req(r),ierr)
      endif
    enddo

    ! Wait for each process to finish sending information
    do k=1,r
      call mpi_wait(req(k),stat,ierr)
    enddo

    ! Adapt time step to ensure stability
    if(Tforce==0)then
      maxtime=max_time
    else
      maxtime=time_step
    endif

    if(Tforce==0)then
      do lid=1,localNodes
        k=localNodesGID(lid)
        id=stackOrder(k)
        rcv=receivers(id)
        th=0.0
        do ks=1,totgrn
          th=th+change_localS(id,ks)
          if(change_localS(id,ks)<0)then
            mtime=-alay_sed(id,ks)/change_localS(id,ks)
            maxtime=min(maxtime,mtime)
          endif
        enddo
        if(tcoordX(id)<minx.or.tcoordX(id)>maxx.or.tcoordY(id)<miny &
            .or.tcoordY(id)>maxy) change_localS(id,1:totgrn)=0.

        if(rcv==id.and.th>0..and.voronoiCell(id)%btype<0)then
          if(watercell(id)==0.)then
            maxh=0.
            do q=1,delaunayVertex(id)%ngbNb
              m=delaunayVertex(id)%ngbID(q)
              maxh=max(maxh,spmZ(m)-spmZ(id))
            enddo
            if(maxh>0.)then
              mtime=maxh/th
              maxtime=min(maxtime,mtime)
            endif
          else
            mtime=watercell(id)/th
            maxtime=min(maxtime,mtime)
          endif
        endif

      enddo
    endif
    call mpi_allreduce(maxtime,dt,1,mpi_double_precision,mpi_min,badlands_world,rc)
    time_step=dt
    if(time_step<force_time) time_step=force_time

    ! Perform elevation and regolith evolution
    do lid=1,localNodes
      id=localNodesGID(lid)
      k=stackOrder(id)
      if(voronoiCell(k)%border==0)then
        if(tcoordX(k)==minx.and.bounds(3)==0)change_local(k)=0.
        if(tcoordX(k)==maxx.and.bounds(4)==0)change_local(k)=0.
        if(tcoordY(k)==miny.and.bounds(2)==0)change_local(k)=0.
        if(tcoordY(k)==maxy.and.bounds(1)==0)change_local(k)=0.
        ! Update surface
        th=0.0
        do ks=1,totgrn
          th=th+time_step*change_localS(k,ks)
          nsed(k,ks)=alay_sed(k,ks)+time_step*change_localS(k,ks)
          if(abs(nsed(k,ks))<1.e-6)then
            if(nsed(k,ks)>0.) th=th-nsed(k,ks)
            nsed(k,ks)=0.
          endif
          if(nsed(k,ks)<-1.e-6)then
            nsed(k,ks)=0.
            th=th-time_step*change_localS(k,ks)-alay_sed(k,ks)
          endif
        enddo
        nZ(k)=spmZ(k)+th
        nth(k)=alay_thick(k)+th
        if(nth(k)<1.e-6)then
          nth(k)=0.
          nsed(k,1:totgrn)=0.
          nZ(k)=spmZ(k)-alay_thick(k)
        endif
      endif
    enddo

  end subroutine geomorphic_evolutionS
  ! =====================================================================================
  subroutine detachmentlimitedS(id,distance,maxh,ST,Qs)

    integer::id,rcv,ks,eroOn
    real(kind=8),dimension(totgrn)::ST,Qs,frac,eroh
    real(kind=8)::dh,distance,maxh,vh


    ! Stream Power Law (detachment-limited) - bedrock incision
    rcv=receivers(id)
    dh=0.95*(spmZ(id)-spmZ(rcv))
    if(dh<0.001)dh=0.
    ST=0.
    eroOn=0
    eroh=0.

    ! Volumic fraction of each sediment class present in the bed
    frac(1:totgrn)=alay_sed(id,1:totgrn)/alay_thick(id)

    ! In case erosion occurs
    if(rcv/=id.and.dh>0..and.watercell(id)<1.e-6.and.spmZ(id)>=gsea%actual_sea)then
      ! If next vertex is potentially erosional
      if(watercell(rcv)<1.e-6)then
        eroOn=1
        vh=min(dh,alay_thick(id))
        do ks=1,totgrn
          eroh(ks)=frac(ks)*vh
          ST(ks)=-frac(ks)*sediments(ks)%ero*(discharge(id))**spl_m*(dh/distance)**spl_n
        enddo
      ! Limit erosion to the elevation difference between current node and lake level
      else
        vh=min(0.99*(spmZ(id)-(watercell(rcv)+spmZ(rcv))),alay_thick(id))
        if(vh>0.)then
          eroOn=1
          do ks=1,totgrn
            eroh(ks)=frac(ks)*vh
            ST(ks)=-frac(ks)*sediments(ks)%ero*(discharge(id))**spl_m*(dh/distance)**spl_n
          enddo
        endif
      endif
    endif

    ! Stability criteria for deposition solutions
    if(eroOn==0.and.perosive==0)then

      ! Marine
      if(maxh>0..and.spmZ(id)<gsea%actual_sea)then

        if(Tforce==0)then
          lp1: do ks=1,totgrn
            if(Qs_inS(id,ks)*max_time/voronoiCell(id)%area<maxh)then
              ST(ks)=Qs_inS(id,ks)/voronoiCell(id)%area
              Qs(ks)=0.
              maxh=maxh-Qs_inS(id,ks)*max_time/voronoiCell(id)%area
            else
              ST(ks)=maxh/max_time
              Qs(ks)=Qs_inS(id,ks)-ST(ks)*voronoiCell(id)%area
              maxh=0.
            endif
            if(maxh<1.e-6) exit lp1
          enddo lp1

        else
          lp2: do ks=1,totgrn
            if(Qs_inS(id,ks)*time_step/voronoiCell(id)%area<maxh)then
              ST(ks)=Qs_inS(id,ks)/voronoiCell(id)%area
              Qs(ks)=0.
              maxh=maxh-Qs_inS(id,ks)*time_step/voronoiCell(id)%area
            else
              ST(ks)=maxh/time_step
              Qs(ks)=Qs_inS(id,ks)-ST(ks)*voronoiCell(id)%area
              maxh=0.
            endif
            if(maxh<1.e-6) exit lp2
          enddo lp2
        endif

      ! Fill depression
      elseif(watercell(id)>0.0001.and.rcv/=id)then
        dh=watercell(id)
        if(Tforce==0)then
          lp3: do ks=1,totgrn
            if(Qs_inS(id,ks)*max_time/voronoiCell(id)%area<dh)then
              ST(ks)=Qs_inS(id,ks)/voronoiCell(id)%area
              Qs(ks)=0.
              dh=dh-Qs_inS(id,ks)*max_time/voronoiCell(id)%area
            else
              ST(ks)=dh/max_time
              Qs(ks)=Qs_inS(id,ks)-ST(ks)*voronoiCell(id)%area
              dh=0.
            endif
            if(dh<1.e-6) exit lp3
          enddo lp3
        else
          lp4: do ks=1,totgrn
            if(Qs_inS(id,ks)*time_step/voronoiCell(id)%area<dh)then
              ST(ks)=Qs_inS(id,ks)/voronoiCell(id)%area
              Qs(ks)=0.
              dh=dh-Qs_inS(id,ks)*time_step/voronoiCell(id)%area
            else
              ST(ks)=dh/time_step
              Qs(ks)=Qs_inS(id,ks)-ST(ks)*voronoiCell(id)%area
              dh=0.
            endif
            if(dh<1.e-6) exit lp4
          enddo lp4
        endif

      elseif(rcv==id)then
        ST(1:totgrn)=Qs_inS(id,1:totgrn)/voronoiCell(id)%area
        Qs=0.

      else
        Qs(1:totgrn)=Qs_inS(id,1:totgrn)
      endif

    ! For purely erosive case
    elseif(eroOn==0.and.perosive==1)then
      Qs(1:totgrn)=Qs_inS(id,1:totgrn)

    ! Stability criteria for erosion solution
    elseif(eroOn==1)then
      if(Tforce==0)then
        do ks=1,totgrn
          if(-ST(ks)*max_time>eroh(ks)) max_time=min(max_time,-eroh(ks)/ST(ks))
          Qs(ks)=-ST(ks)*voronoiCell(id)%area+Qs_inS(id,ks)
        enddo
      else
        do ks=1,totgrn
          if(-ST(ks)*time_step>eroh(ks)) ST(ks)=-eroh(ks)/time_step
          Qs(ks)=-ST(ks)*voronoiCell(id)%area+Qs_inS(id,ks)
        enddo
      endif
    endif

  end subroutine detachmentlimitedS
  ! =====================================================================================
end module stratmorph
