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
!       Filename:  StratalFunctions.f90
!
!    Description:  Stratigraphic grid evolution functions.
!
!        Version:  1.0
!        Created:  03/06/15 09:05:27
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module strata_evol

  use parallel
  use parameters
  use out_stratal
  use stratal_read
  use stratal_class
  use external_forces

  implicit none

contains

  ! =====================================================================================
  subroutine StrataGen

    ! Get stratal definition
    call stratal_parser

    ! Built stratigraphy grid
    if(totgrn>0) call constructStrata

    return

  end subroutine StrataGen
  ! =====================================================================================
  subroutine compute_stratal_displacement

    integer::k,id

    do k=1,upartN
      id=unodeID(k)
      lay_base(k)=lay_base(k)+tvertDisp(id)*time_step
    enddo

  end subroutine compute_stratal_displacement
  ! =====================================================================================
  subroutine update_stratigraphy_layer

    integer::k,gid,ks,ll

    real(kind=8)::remain,sed

    ! Record the active layer changes in the stratigraphic mesh
    do k=1,upartN
      gid=unodeID(k)
      do ks=1,totgrn
        sed=alay_sed(gid,ks)-alay_dsed(gid,ks)
        ! Deposition
        if(sed>=1.e-6)then
          lay_sed(k,layNb,ks)=lay_sed(k,layNb,ks)+sed
          lay_thick(k,layNb)=lay_thick(k,layNb)+sed
        elseif(sed<-1.e-6)then
          remain=-sed
          lpp: do ll=layNb,1,-1
            if(lay_sed(k,ll,ks)>remain)then
              lay_sed(k,ll,ks)=lay_sed(k,ll,ks)-remain
              lay_thick(k,ll)=lay_thick(k,ll)-remain
              exit lpp
            else
              remain=remain-lay_sed(k,ll,ks)
              lay_thick(k,ll)=lay_thick(k,ll)-lay_sed(k,ll,ks)
              lay_sed(k,ll,ks)=0.0
            endif
            if(lay_thick(k,ll)<1.e-6)then
              lay_thick(k,ll)=0.0
              lay_sed(k,ll,1:totgrn)=0.0
            endif
          enddo lpp
          if(ll>=1)then
            if(lay_thick(k,ll)<1.e-6)then
              lay_thick(k,ll)=0.0
              lay_sed(k,ll,1:totgrn)=0.0
            endif
          endif
        endif
      enddo
    enddo

    ! Build new active layer
    if(simulation_time<time_end) call buildActiveLayer

    return

  end subroutine update_stratigraphy_layer
  ! =====================================================================================
  subroutine buildActiveLayer

    integer::k,l,ks,gid

    real(kind=8)::th,sed(totgrn),prop(totgrn)
    real(kind=8),dimension(dnodes)::tsed,gsed

    if(.not.allocated(alay_thick))then
      allocate(alay_thick(dnodes))
      allocate(alay_sed(dnodes,totgrn))
      allocate(alay_dsed(dnodes,totgrn))
    endif

    ! Get the composition and thickness of the active layer
    alay_thick=0.0
    alay_sed=0.0
    do k=1,upartN
      th=0.0
      sed=0.0
      gid=unodeID(k)
      lp:do l=layNb,1,-1
        th=th+lay_thick(k,l)
        if(th>active_thick)then
          th=active_thick
          do ks=1,totgrn
            prop(ks)=lay_sed(k,l,ks)/lay_thick(k,l)
            sed(ks)=sed(ks)+prop(ks)*(th-alay_thick(gid))
          enddo
          alay_thick(gid)=th
        else
          alay_thick(gid)=th
          do ks=1,totgrn
            sed(ks)=sed(ks)+lay_sed(k,l,ks)
          enddo
        endif
        if(th>=active_thick)then
          alay_sed(gid,1:totgrn)=sed(1:totgrn)
          exit lp
        endif
      enddo lp
    enddo

    ! Broadcast local stratigraphic layer globally
    call mpi_allreduce(mpi_in_place,alay_thick,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
    do ks=1,totgrn
      tsed(1:dnodes)=alay_sed(1:dnodes,ks)
      call mpi_allreduce(tsed,gsed,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
      alay_sed(1:dnodes,ks)=gsed(1:dnodes)
      alay_dsed(1:dnodes,ks)=gsed(1:dnodes)
    enddo

    return

  end subroutine buildActiveLayer
  ! =====================================================================================
end module strata_evol
