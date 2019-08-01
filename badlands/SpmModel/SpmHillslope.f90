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
!       Filename:  SpmHillslope.f90
!
!    Description:  Defines the hillslope transport mechanism by diffusion
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module hillslope

  use parallel
  use topology
  use parameters
  use hydroUtil
  use hydrology
  use external_forces

  implicit none

contains

  ! =====================================================================================
  subroutine CFLdiffusion

    integer::lid,k,n,p,id
    real(kind=8)::distance,denom,Cdiff,ldt(2),dt(2)

  	time_step=CFL_diffusion

    do lid=1,localNodes
      n=localNodesGID(lid)
      k=stackOrder(n)
      if(voronoiCell(k)%border==0)then
      	do p=1,delaunayVertex(k)%ngbNb
      		if(delaunayVertex(k)%ngbID(p)>0)then
            if(voronoiCell(delaunayVertex(k)%ngbID(p))%border==0)then
        			id=delaunayVertex(k)%ngbID(p)
        			distance=delaunayVertex(k)%distance(p)
        			! CFL linear diffusion
              Cdiff=max(Cdiffusion(1),Cdiffusion(2))
              if(simulation_time==time_start.and.Cdiff>0.)then
        				CFL_diffusion=min(CFL_diffusion,0.1*(distance)**2/(2.0*Cdiff))
        				time_step=min(time_step,CFL_diffusion)
        			endif
              if(simulation_time>time_start.and.update3d.and.Cdiff>0.)then
                CFL_diffusion=min(CFL_diffusion,0.1*(distance)**2/(2.0*Cdiff))
                time_step=min(time_step,CFL_diffusion)
              endif
        			! CFL depth-dependent diffusion
              Cdiff=max(Cdiffusion_d(1),Cdiffusion_d(2))
        			if(Cdiff>0.)then
        				denom=Cdiffusion_d(1)*spmH(k)**Cdiff_m*(abs(spmZ(k)-spmZ(id))/distance)**(Cdiff_n-1.0)
                time_step=min(time_step,distance/denom)
        			endif
        			! CFL non-linear diffusion
              Cdiff=max(Cdiffusion_nl(1),Cdiffusion_nl(2))
        			if(Cdiff>0.)then
                time_step=min(time_step,0.1*(distance**2-abs(spmZ(k)-spmZ(id))/slope_critical**2) &
                  /(2.0*Cdiff))
              endif
        		endif
          endif
      	enddo
      endif
    enddo

    if(time_step>display_interval) time_step=display_interval
    ldt(1)=CFL_diffusion
    ldt(2)=time_step
    call mpi_allreduce(ldt,dt,2,mpi_double_precision,mpi_min,badlands_world,rc)
    CFL_diffusion=dt(1)
    time_step=dt(2)

  end subroutine CFLdiffusion
  ! =====================================================================================
  subroutine hillslope_flux(id,LDL,DDD,NDL,diffH)

    integer::id,n,p
    real(kind=8)::LDL,DDD,NDL,dh,edge,distance,diffH

    LDL=0.
    DDD=0.
    NDL=0.
    diffH=1.e6
    do p=1,delaunayVertex(id)%ngbNb
      if(delaunayVertex(id)%ngbID(p)>0)then
        if(voronoiCell(delaunayVertex(id)%ngbID(p))%border==0)then
      		n=delaunayVertex(id)%ngbID(p)
      		! Define neighborhood parameters
      		dh=(spmZ(n)-spmZ(id))
          if(dh>=0.) diffH=min(dh,diffH)
      		edge=delaunayVertex(id)%voronoi_edge(p)
      		distance=delaunayVertex(id)%distance(p)
      		! Hillslope transport by linear diffusion
      		if(Cdiffusion(1)>0..or.Cdiffusion(2)>0.) LDL=LDL+ &
            edge*dh/distance
      		! Hillslope transport by depth-dependent diffusion
      		if(Cdiffusion_d(1)>0..or.Cdiffusion_d(2)>0.) DDD=DDD+ &
            edge*(dh/distance)**Cdiff_n
      		! Hillslope transport by non-linear diffusion
      		if(Cdiffusion_nl(1)>0..or.Cdiffusion_nl(2)>0.) NDL=NDL+ &
            (edge*dh/distance)/(1.0-(abs(dh)*slope_critical/distance))
      	endif
      endif
    enddo
    if(diffH>9.99e5) diffH=0.
    if(spmZ(id)<gsea%actual_sea)then
      diffH=gsea%actual_sea-spmZ(id)
    elseif(watercell(id)>0.)then
      diffH=min(watercell(id),diffH)
    endif
    diffH=0.95*diffH

    ! Update elevation change for each hillslope diffusion-like process
    if(spmZ(id)<gsea%actual_sea.or.watercell(id)>0.)then
      LDL=LDL*Cdiffusion(2)/voronoiCell(id)%area
      DDD=DDD*Cdiffusion_d(2)*(spmH(id))**Cdiff_m/voronoiCell(id)%area
      NDL=NDL*Cdiffusion_nl(2)/voronoiCell(id)%area
    else
      LDL=LDL*Cdiffusion(1)/voronoiCell(id)%area
      DDD=DDD*Cdiffusion_d(1)*(spmH(id))**Cdiff_m/voronoiCell(id)%area
      NDL=NDL*Cdiffusion_nl(1)/voronoiCell(id)%area
    endif

  end subroutine hillslope_flux
  ! =====================================================================================
end module hillslope
