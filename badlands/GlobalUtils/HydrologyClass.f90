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
!       Filename:  HydrologyClass.f90
!
!    Description:  Store hydrological parameters
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module hydroUtil

  use parallel

  implicit none

  ! Logical flag for regolith file
  logical::regofileflag

  ! Receivers nodes
  integer,dimension(:),allocatable::receivers

  ! Number of donors
  integer,dimension(:),allocatable::donorsList

  ! Index array
  integer,dimension(:),allocatable::indexArray

  ! Donors information array
  integer,dimension(:,:),allocatable::donorInfo

  ! Nodes ordering stack
  integer,dimension(:),allocatable::stackOrder,lstackOrder

  ! List of baselevel nodes
  integer::baseNb
  integer,dimension(:),allocatable::baselist

  ! Discharge
  real(kind=8),dimension(:),allocatable::discharge

  ! Chi parameter
  real(kind=8),dimension(:),allocatable::chi

  ! Chi parameter
  integer,dimension(:),allocatable::bsID

  ! Maximum number of receiver nodes
  integer::maxrcvs

  ! Stream network accumulation minimal value
  real(kind=8)::accu_thres

  ! Strahler stream order
  integer,dimension(:),allocatable::strahler

  ! Cathcment ID
  integer::catchID
!   integer,dimension(:),allocatable::catchmentID
  integer,dimension(:),allocatable::subcatchmentID
  integer,dimension(:),allocatable::subcatchmentProc

  ! Simulation time parameters
  integer::Tforce
  real(kind=8)::time_step,simulation_time,udw_time,force_time
  real(kind=8)::time_start,time_end,time_display,cpl1_time,cpl2_time,cpl3_time,cpl4_time,cpl5_time
  real(kind=8)::display_interval,CFL_diffusion

  ! Diffusion coefficient
  real(kind=8),dimension(2)::Cdiffusion ! linear
  real(kind=8),dimension(2)::Cdiffusion_d ! depth-dependent
  real(kind=8)::Cdiff_m ! depth-dependent
  real(kind=8)::Cdiff_n ! depth-dependent
  real(kind=8),dimension(2)::Cdiffusion_nl ! non-linear
  real(kind=8)::slope_critical

  ! Surface erodibility coefficient
  integer::perosive
  real(kind=8)::Cerodibility

  ! Ice erodibility coefficient
  real(kind=8)::IceEro

  ! Stream Power Law coefficient
  real(kind=8)::spl_m
  real(kind=8)::spl_n

  ! Sefiment transport efficiency coefficient
  real(kind=8)::Cefficiency

  ! Fraction of total load delivered to channels as bedload
  real(kind=8)::Fracbed

  ! Stream Power Law coefficient
  real(kind=8)::stl_m
  real(kind=8)::stl_n

  ! Regolith formation
  real(kind=8)::regoProd
  real(kind=8)::regoDepth
  real(kind=8)::soil_density
  real(kind=8)::rock_density

  ! Capacity model coefficient
  real(kind=8)::chan_exp
  real(kind=8)::chan_width
  real(kind=8)::bed_length
  real(kind=8)::sed_length
  real(kind=8)::stream_ero
  real(kind=8)::bed_sed_interface

  ! New elevation change
  real(kind=8),dimension(:),allocatable::newZ,spmZ,spmH

  ! Sediment influx
  real(kind=8),dimension(:),allocatable::Qs_in

  ! Subcatchement declaration type
  integer,dimension(:),allocatable::subcatchNb

  ! Network parallelisation stream declaration
  integer::localNodes
  integer,dimension(:),allocatable::sendprocID,rcvprocNb
  integer,dimension(:,:),allocatable::rcvsendID,rcvprocID
  integer,dimension(:),allocatable::localNodesGID

end module hydroUtil
