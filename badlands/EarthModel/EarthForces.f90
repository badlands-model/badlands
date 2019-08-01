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
!       Filename:  EarthForces.f90
!
!    Description:  Defines the geodynamic evolution and the rainfall.
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module earthforces

  use parallel
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none

  integer,dimension(:),allocatable::surfID

contains

  ! =====================================================================================
  subroutine displacement

    integer::kn,k,iu,p,id,m,n


    if(simulation_time==time_start.and..not.udwFlag)then
      if(.not.allocated(disp_fill)) allocate(disp_fill(disp%event))
      if(disp_time(1,1)>time_start)vdisp%event=1
      ! Find the vertical displacement event number
      do k=1,disp%event
        vdisp%event=vdisp%event+1
        if(disp_time(k,1)>=disp_time(k,2))then
          if(pet_id==0)print*,'ERROR: reading displacements time declaration in event:',k
          call mpi_finalize(rc)
        endif
        if(k<disp%event)then
          if(disp_time(k,2)>disp_time(k+1,1))then
            if(pet_id==0)print*,'ERROR: reading displacements time declaration between events:',k,k+1
            call mpi_finalize(rc)
          endif
          if(disp_time(k,2)/=disp_time(k+1,1)) vdisp%event=vdisp%event+1
        endif
      enddo

      ! If last displacement period stops before the end of the simulation
      ! add a new event during the last period
      if(disp_time(disp%event,2)<time_end) vdisp%event=vdisp%event+1

      ! In case the number of vertical displacement is different from the number of
      ! of tectonic events defined in the XmL input file defines the new event parameters
      if(vdisp%event/=disp%event)then
        allocate(vdisp_time(vdisp%event,2),vdisp_fill(vdisp%event))
        kn=1
        if(disp_time(1,1)>time_start)then
          vdisp_time(kn,1)=time_start
          vdisp_time(kn,2)=disp_time(1,1)
          vdisp_fill(kn)=0
          kn=kn+1
        endif

        do k=1,disp%event
          vdisp_time(kn,1)=disp_time(k,1)
          vdisp_time(kn,2)=disp_time(k,2)
          vdisp_fill(kn)=k
          kn=kn+1
          if(k<disp%event)then
            if(disp_time(k,2)/=disp_time(k+1,1))then
              vdisp_time(kn,1)=disp_time(k,2)
              vdisp_time(kn,2)=disp_time(k+1,1)
              vdisp_fill(kn)=0
              kn=kn+1
            endif
          endif
        enddo

        if(disp_time(disp%event,2)<time_end)then
          vdisp_time(kn,2)=time_end
          vdisp_time(kn,1)=disp_time(disp%event,2)
          vdisp_fill(kn)=0
        endif

        deallocate(disp_time,disp_fill)
        allocate(disp_time(vdisp%event,2),disp_fill(vdisp%event))
        disp%event=vdisp%event
        do k=1,disp%event
          disp_time(k,1)=vdisp_time(k,1)
          disp_time(k,2)=vdisp_time(k,2)
          disp_fill(k)=vdisp_fill(k)
        enddo
        deallocate(vdisp_time,vdisp_fill)

      ! Otherwise define the displacement event number directly
      else
        do k=1,disp%event
          disp_fill(k)=k
        enddo
      endif

    elseif(simulation_time==time_start.and.udwFlag)then
      if(.not.allocated(disp_fill)) allocate(disp_fill(disp%event))
      do k=1,disp%event
        disp_fill(k)=k
      enddo
    endif

    ! Find actual event number
    disp%actual=0
    do k=1,disp%event
      if(simulation_time<disp_time(k,2).and.simulation_time>=disp_time(k,1)) disp%actual=k
    enddo

    ! Open displacements field file for the considered event number
    kn=disp%actual
    if(disp_fill(kn)>0)then
      iu=54
      open(iu,file=fgeodyn(disp_fill(kn)),status="old",action="read",iostat=rc)
      if(rc/=0)then
        print*,'Failed to open tectonic file'
        call mpi_finalize(rc)
      endif

      if(udwFlag)then
        if(.not.allocated(surfID))then
          allocate(surfID(nx*ny))
          id=0
          m=0
          surfID=0
          do k=1,ny+2
            do p=1,nx+2
              m=m+1
              if(p>1.and.p<nx+2.and.k>1.and.k<ny+2)then
                id=id+1
                surfID(id)=m
              endif
            enddo
          enddo
        endif
      endif

      m=0
      do k=1,ny+2
        do p=1,nx+2
          m=m+1
          if(p>1.and.p<nx+2.and.k>1.and.k<ny+2)then
            if(udwFlag)then
              if(disp3d)then
                read(iu,*)n,rDisp(surfID(n),1:3)
              else
                read(iu,*)n,rDisp(surfID(n),1)
              endif
            else
              if(disp3d)then
                read(iu,*)rDisp(m,1:3)
              else
                read(iu,*)rDisp(m,1)
              endif
            endif
          endif
        enddo
      enddo
      close(iu)
      ! Get the regular grid corners values
      rDisp(1,:)=rDisp(nx+4,:)
      rDisp(nx+2,:)=rDisp(2*nx+3,:)
      rDisp(bnbnodes-nx-1,:)=rDisp(bnbnodes-2*(nx+2)+2,:)
      rDisp(bnbnodes,:)=rDisp(bnbnodes-(nx+2)-1,:)
      ! Get the South/North displacement values
      id=0
      do k=1,ny+2
        do p=1,nx+2
          id=id+1
          if(k==1.and.p>1.and.p<nx+2)then
            rDisp(id,:)=rDisp(id+nx+2,:)
          elseif(k==ny+2.and.p>1.and.p<nx+2)then
            rDisp(id,:)=rDisp(id-nx-2,:)
          endif
        enddo
      enddo
      ! Get the East/West displacement values
      id=0
      do k=1,ny+2
        do p=1,nx+2
          id=id+1
          if(p==1.and.k>1.and.k<ny+2)then
            rDisp(id,:)=rDisp(id+1,:)
          elseif(p==nx+2.and.k>1.and.k<ny+2)then
            rDisp(id,:)=rDisp(id-1,:)
          endif
        enddo
      enddo
    ! If this event doesn't correspond to any defined one initialise displacement values to 0
    else
      rDisp=0.0
    endif

    ! Update displacement rate
    do k=1,bnbnodes
      if(disp3d)then
        if(simulation_time==time_start.and.disp_time(disp%actual,1)<time_start)then
          rhxDisp(k)=0.0
          rhyDisp(k)=0.0
          rvertDisp(k)=0.0
        else
          if(disp%disptime==0.0)then
            rhxDisp(k)=rDisp(k,1) !/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
            rhyDisp(k)=rDisp(k,2) !/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
            rvertDisp(k)=rDisp(k,3) !/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
          else
            rhxDisp(k)=rDisp(k,1)*disp%disptime/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
            rhyDisp(k)=rDisp(k,2)*disp%disptime/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
            rvertDisp(k)=rDisp(k,3)*disp%disptime/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
          endif
        endif
      else
        rvertDisp(k)=rDisp(k,1)/(disp_time(disp%actual,2)-disp_time(disp%actual,1))
      endif
    enddo

    ! Define next coupling time
    if(.not.disp3d)then
      if(disp%actual<disp%event)then
        cpl2_time=disp_time(disp%actual+1,1)
      else
        cpl2_time=time_end+1000.
      endif
    elseif(disp3d.and.disp%disptime>0..and.disp%actual>0)then
      if(disp%actual==disp%event.and.simulation_time<disp_time(disp%actual,2))then
        cpl2_time=cpl2_time+disp%disptime
      elseif(disp_time(disp%actual+1,1)>cpl2_time+disp%disptime)then
        cpl2_time=cpl2_time+disp%disptime
      else
        cpl2_time=disp_time(disp%actual+1,1)
      endif
    endif

  end subroutine displacement
  ! =====================================================================================
  subroutine rainfall

    integer::n,event_nb,iu,p,k,id,m

    ! Find the rain event if it exists and update it if required
    event_nb=0
    find_event: do n=1,rain_event
      if(rain_tstart(n)<=simulation_time.and.rain_tend(n)>simulation_time)then
        event_nb=n
        exit find_event
      endif
    enddo find_event
    ! Define next coupling time for rain map
    if(event_nb<rain_event)then
      cpl1_time=rain_tstart(event_nb+1)
    else
      cpl1_time=time_end+1000.
    endif
    ! Is there a new map to load
    if(.not.allocated(rainVal)) allocate(rainVal(bnbnodes))
    rainVal=0.
    if(event_nb>0)then
      iu=51
      if(frainval(event_nb)<0.)then
        open(iu,file=frainmap(event_nb),status="old",action="read",iostat=rc)
        if(rc/=0)then
          print*,'Failed to open namelist file rain map file'
          call mpi_finalize(rc)
        endif
        id=0
        m=0
        do k=1,ny+2
          do p=1,nx+2
            m=m+1
            if(p>1.and.p<nx+2.and.k>1.and.k<ny+2)then
              id=id+1
              read(iu,*)rainVal(m)
            endif
          enddo
        enddo
        close(iu)
        ! Get the regular grid corners values
        rainVal(1)=rainVal(nx+4)
        rainVal(nx+2)=rainVal(2*nx+3)
        rainVal(bnbnodes-nx-1)=rainVal(bnbnodes-2*(nx+2)+2)
        rainVal(bnbnodes)=rainVal(bnbnodes-(nx+2)-1)
        ! Get the South/North rain values
        id=0
        do k=1,ny+2
          do p=1,nx+2
            id=id+1
            if(k==1.and.p>1.and.p<nx+2)then
              rainVal(id)=rainVal(id+nx+2)
            elseif(k==ny+2.and.p>1.and.p<nx+2)then
              rainVal(id)=rainVal(id-nx-2)
            endif
          enddo
        enddo
        ! Get the East/West rain values
        id=0
        do k=1,ny+2
          do p=1,nx+2
            id=id+1
            if(p==1.and.k>1.and.k<ny+2)then
              rainVal(id)=rainVal(id+1)
            elseif(p==nx+2.and.k>1.and.k<ny+2)then
              rainVal(id)=rainVal(id-1)
            endif
          enddo
        enddo
      else
        rainVal=frainval(event_nb)
      endif
    endif

  end subroutine rainfall
  ! =====================================================================================
end module earthforces
