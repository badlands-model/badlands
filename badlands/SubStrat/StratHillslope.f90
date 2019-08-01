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
!       Filename:  StratHillslope.f90
!
!    Description:  Defines the hillslope transport mechanism by diffusion for stratigraphic grid
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module strat_hillslope

  use parallel
  use topology
  use parameters
  use hydroUtil
  use hydrology
  use stratal_class
  use external_forces

  implicit none

contains

  ! =====================================================================================
  subroutine CFLdiffusionS

    integer::lid,k,n,p,id
    real(kind=8)::distance,Cdiff,ldt(2),dt(2)

  	time_step=CFL_diffusion

    Cdiff=0.0
    do k=1,totgrn
      Cdiff=max(Cdiff,sediments(k)%diffa)
      Cdiff=max(Cdiff,sediments(k)%diffm)
    enddo

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
              if(simulation_time==time_start.and.Cdiff>0.)then
        				CFL_diffusion=min(CFL_diffusion,0.1*(distance)**2/(2.0*Cdiff))
        				time_step=min(time_step,CFL_diffusion)
        			endif
              if(simulation_time>time_start.and.update3d.and.Cdiff>0.)then
                CFL_diffusion=min(CFL_diffusion,0.1*(distance)**2/(2.0*Cdiff))
                time_step=min(time_step,CFL_diffusion)
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

  end subroutine CFLdiffusionS
  ! =====================================================================================
  subroutine hillslope_fluxS(id,LDL,diffH)

    integer::id,n,p,ks

    real(kind=8),dimension(totgrn)::LDL,frac
    real(kind=8)::dh,edge,distance,diffH

    LDL=0.
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
          if(dh>=0.and.alay_thick(n)>0.0)then
            frac(1:totgrn)=alay_sed(n,1:totgrn)/alay_thick(n)
          elseif(alay_thick(id)>0.0)then
            frac(1:totgrn)=alay_sed(id,1:totgrn)/alay_thick(id)
          elseif(dh>=0..and.alay_thick(n)==0.0)then
            frac(1:totgrn)=alay_sed(id,1:totgrn)/alay_thick(id)
          elseif(alay_thick(id)==0.0)then
            print*,'Issue in multi-size sediment diffusion.',id,n
          endif
          do ks=1,totgrn
        		LDL(ks)=LDL(ks)+edge*dh*frac(ks)/distance
          enddo
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
      do ks=1,totgrn
        LDL(ks)=LDL(ks)*sediments(ks)%diffm/voronoiCell(id)%area
      enddo
    else
      do ks=1,totgrn
        LDL(ks)=LDL(ks)*sediments(ks)%diffa/voronoiCell(id)%area
      enddo
    endif

  end subroutine hillslope_fluxS
  ! =====================================================================================
end module strat_hillslope
