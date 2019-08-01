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
!       Filename:  WaveModel.f90
!
!    Description:  SWAN wave model.
!
!        Version:  1.0
!        Created:  06/10/15 20:02:47
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module wave_model

  use parallel
  use topology
  use parameters
  use external_forces
  use ocean_circ
  use swan_coupler_export
  use swan_coupler_parallel
  use swan_coupler_functions

  implicit none

  real,dimension(6)::hcast
  type(badlands_WADecomp),dimension(:),allocatable::wdcp

contains

  ! =====================================================================================
  subroutine wave_initialisation

    if(waveON==0) return

    call swan_createinput
    call swan_initialize(badlands_world,swaninput(1:80),swaninfo(1:80))
    if(allocated(wdcp)) deallocate(wdcp)
    allocate(wdcp(npets))
    call swan_decomposition_grid(npets,wdcp)

    if(allocated(waveU)) deallocate(waveU)
    allocate(waveU(nx+2,ny+2,season))
    if(allocated(waveV)) deallocate(waveV)
    allocate(waveV(nx+2,ny+2,season))

    return

  end subroutine wave_initialisation
  ! ============================================================================
  subroutine wave_update_bathymetry

    integer::p,ic,jc,ii,jj,ks,i,j

    real,dimension(nx,ny)::tempZ

    if(waveON==0) return

    ! Read changes in topographic regular grid.
    p=0
    jc=1
    ic=0
    tempZ=0.
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          ic=ic+1
          tempZ(ic,jc)=real(rcoordZ(p))
          if(i==nx+1)then
            ic=0
            jc=jc+1
          endif
        endif
      enddo
    enddo

    ! Get bathymetry from badlands
    ii=0
    ks=pet_id+1
    bathyfield=0.
    do i=wdcp(ks)%minX_id,wdcp(ks)%maxX_id
      ii=ii+1
      jj=0
      do j=wdcp(ks)%minY_id,wdcp(ks)%maxY_id
        jj=jj+1
        bathyfield(ii,jj)=real(gsea%actual_sea)-tempZ(i,j)
        if(bathyfield(ii,jj)<0.) bathyfield(ii,jj)=-999999.0
        if(bathyfield(ii,jj)>wavebase) bathyfield(ii,jj)=-999999.0
      enddo
    enddo

    return

  end subroutine wave_update_bathymetry
  ! ============================================================================
  subroutine swan_createinput

    integer::ios,iu
    real::mx_x,mx_y
    character(len=128)::stg1,stg2

    ! Create the input for SWAN simulation
    swaninput='swan.swn'
    call noblnk(swaninput)
    call addpath1(swaninput)
    swaninfo='swanInfo'
    call noblnk(swaninfo)
    call addpath1(swaninfo)
    swanbot='swan.bot'
    call noblnk(swanbot)
    call addpath1(swanbot)
    swanout='swan.csv'
    call noblnk(swanout)
    call addpath2(swanout)
    if(pet_id==0)then
      iu=342
      mx_x=real(dx)*(nx-1)
      mx_y=real(dx)*(ny-1)
      open(iu,file=swaninput,status="replace",action="write",iostat=ios)
      write(iu,'(a17)') "PROJECT ' ' 'S01' "
      write(iu,102) "CGRID REG",minx,miny,0,mx_x,mx_y,nx-1,ny-1,"CIRCLE 36 0.05 1.0"
102   format(a9,1x,f12.3,1x,f12.3,1x,i1,1x,f12.3,1x,f12.3,1x,i4,1x,i4,1x,a19)
      write(iu,103) "INPGRID BOTTOM REG",minx,miny,0,1,1,mx_x,mx_y,'EXC -999999.000'
103   format(a18,1x,f12.3,1x,f12.3,1x,i1,1x,i1,1x,i1,1x,f12.3,1x,f12.3,1x,a15)
      stg1="READINP BOTTOM 1 '"
      call append_str2(stg1,swanbot)
      stg2="' 3 0 FREE"
      call append_str2(stg1,stg2)
      write(iu,*)trim(stg1)
      write(iu,'(a42)') "BOUNd SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGR"
      write(iu,'(a7)') 'DIFFRAC'
      write(iu,'(a8)') 'FRICTION'
      ! write(iu,'(a26)') 'BREAKING CONSTANT 1.1 0.73'
      ! write(iu,'(a8)') 'WCAPPING'
      ! write(iu,'(a8)') 'OFF QUAD'
      write(iu,'(a4)') 'GEN1'
      ! write(iu,'(a4)') 'QUAD'
      ! write(iu,'(a10)') 'GEN3 AGROW'
      ! write(iu,'(a10)') 'OFF BNDCHK'
      write(iu,'(a15,1x,i1,1x,i4,1x,i1,1x,i4)')"GROUP 'gf' SUBG",0,nx-1,0,ny-1
      stg1="TABLE 'gf' IND '"
      call append_str2(stg1,swanout)
      stg2="' XP YP DIR UBOT" !HS PER WLEN UBOT"
      call append_str2(stg1,stg2)
      write(iu,*)trim(stg1)
      write(iu,'(a5,2f12.3)') 'WIND ',forecast_param(5:6)
      ! write(iu,'(a9,4f12.3)')'INIT PAR ',forecast_param(1:4)
      write(iu,'(a7)') 'COMPUTE'
      close(iu)
      open(iu,file=swanbot,status="replace",action="write",iostat=ios)
      write(iu,'(a6)')'-1 -1'
      write(iu,'(a6)')'-1 -1'
      close(iu)
    endif
    call mpi_barrier(badlands_world,rc)

    return

    end subroutine swan_createinput
    ! ============================================================================
    subroutine waves_circulation_run(sgp)

      integer::sgp

      if(waveON==0) return

      ! Define forcing waves parameters
      call import_bathymetry
      hcast(1:6)=forecast_param(1:6)
      call swan_run(hcast)
      call get_wave_velocity(sgp)
      return

    end subroutine waves_circulation_run
    ! ============================================================================
    subroutine get_wave_velocity(sgp)

      integer::i,j,xo,yo,x1,y1,xpart,ypart,ii,jj,p,sgp
      real::pp
      real,dimension(nx,ny)::velwaveV,velwaveU
      real,dimension(nx*ny)::tempU,tempV,gtempU,gtempV

      xo=wdcp(pet_id+1)%minX_id
      yo=wdcp(pet_id+1)%minY_id
      x1=wdcp(pet_id+1)%maxX_id
      y1=wdcp(pet_id+1)%maxY_id
      xpart=wdcp(pet_id+1)%X_nb
      ypart=wdcp(pet_id+1)%Y_nb
      pp=4.*atan(1.)

      do i=1,xpart
        do j=1,ypart
          ! Wave orbital velocity
          if(bathyfield(i,j)<0) exportSWANfields(i,j,1)=0.0
          if(bathyfield(i,j)<0) exportSWANfields(i,j,2)=0.0
        enddo
      enddo

      ii=0
      velwaveU=0.
      velwaveV=0.
      do i=xo,x1
        ii=ii+1
        jj=0
        do j=yo,y1
          jj=jj+1
          velwaveU(i,j)=exportSWANfields(ii,jj,1)
          velwaveV(i,j)=exportSWANfields(ii,jj,2)
        enddo
      enddo

      tempU=-1.e8
      tempV=-1000.
      p=0
      do i=1,nx
        do j=1,ny
          p=p+1
          tempU(p)=velwaveU(i,j)
          tempV(p)=velwaveV(i,j)
        enddo
      enddo

      call mpi_allreduce(tempU,gtempU,nx*ny,real_type,max_type,badlands_world,rc)
      call mpi_allreduce(tempV,gtempV,nx*ny,real_type,max_type,badlands_world,rc)

      p=0
      do i=1,nx+2
        do j=1,ny+2
          if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
            p=p+1
            if(gtempU(p)>0.0)then
              waveU(i,j,sgp)=real(gtempU(p)*cos(gtempV(p)*pp/180.))
              waveV(i,j,sgp)=real(gtempU(p)*sin(gtempV(p)*pp/180.))
            else
              waveU(i,j,sgp)=0.
              waveV(i,j,sgp)=0.
            endif
          endif
        enddo
      enddo

      ! Update border
      waveU(2:nx+1,1,sgp)=waveU(2:nx+1,2,sgp)
      waveU(2:nx+1,ny+2,sgp)=waveU(2:nx+1,ny+1,sgp)
      waveU(1,2:ny+1,sgp)=waveU(2,2:ny+1,sgp)
      waveU(nx+2,2:ny+1,sgp)=waveU(nx+1,2:ny+1,sgp)
      waveV(2:nx+1,1,sgp)=waveV(2:nx+1,2,sgp)
      waveV(2:nx+1,ny+2,sgp)=waveV(2:nx+1,ny+1,sgp)
      waveV(1,2:ny+1,sgp)=waveV(2,2:ny+1,sgp)
      waveV(nx+2,2:ny+1,sgp)=waveV(nx+1,2:ny+1,sgp)

      ! Update corner
      waveU(1,1,sgp)=waveU(2,2,sgp)
      waveU(1,ny+2,sgp)=waveU(2,ny+1,sgp)
      waveU(nx+2,1,sgp)=waveU(nx+1,2,sgp)
      waveU(nx+2,ny+2,sgp)=waveU(nx+1,ny+1,sgp)
      waveV(1,1,sgp)=waveV(2,2,sgp)
      waveV(1,ny+2,sgp)=waveV(2,ny+1,sgp)
      waveV(nx+2,1,sgp)=waveV(nx+1,2,sgp)
      waveV(nx+2,ny+2,sgp)=waveV(nx+1,ny+1,sgp)

      return

    end subroutine get_wave_velocity
    ! =====================================================================================
end module wave_model
