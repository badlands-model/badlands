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
!       Filename:  GlobalClass.f90
!
!    Description:  Store global application parameters
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module parallel

#include "mpif.h"

  integer,parameter::badlands_world=mpi_comm_world

end module parallel

module parameters

  use parallel

  implicit none

  ! MPI type communicator
  integer::real_type,max_type

  logical::updateSPM_elevation,oceanFlag,geodynamicFlag,restartFlag,udwFlag,update3d

  ! Step for restarting simulation
  integer::restartStep,restartPet

  ! Persistent Execution Threads (ID and total number)
  integer::pet_id,npets,rc

  ! Grid / Coupler Component Names
!   character(len=128)::ocean,spm,earth

  ! Number of total grains
  integer::totgrn

  ! Simulation directories and file names
  character(len=128)::regularfile,outdir,xmlfile
  character(len=128)::outdir1,outdir2,runfiles,outputs

  ! Underworld files and folders
  character(len=128)::outdir3,fudw,fudisp,maestro

  ! Regolith file name
  character(len=128)::regofile

  ! Regular structured grid nodes on X/Y directions
  integer::nx,ny

  ! Regular structured grid resolution
  real(kind=8)::dx

  ! Regular structured grid region extensions
  real(kind=8)::minx,miny,maxx,maxy

  ! Regular structured grid X,Y,Z coordinates
  real(kind=8),dimension(:),allocatable::coordX,coordY,coordZ

  ! Regular structured grid X,Y,Z coordinates with added cells
  ! on the edges and vertical displacement
  real(kind=8),dimension(:),allocatable::rcoordX,rcoordY,rcoordZ

  ! Regular structured grid bilinear
  real,dimension(:),allocatable::bilinearX,bilinearY
  real,dimension(:,:),allocatable::bilinearV,bilinearHx,bilinearHy

  ! Conforming Delaunay Triangulation (CDT) parameters
  integer::tnodes,dnodes,delem,dedge

  ! Output parameters for CDT and drainage
  integer::delemo,delemoo,drainOde
  integer,dimension(:),allocatable::outelem,outnode,dglbID,doutelem,doutnode

  ! Voronoi Diagram of the CDT parameters
  integer::vnodes,vedge,vcellIN,velemIN

  ! Conforming Delaunay Triangulation X,Y,Z coordinates and vertical displacement
  real(kind=8),dimension(:),allocatable::tcoordX,tcoordY,tcoordZ

  ! Voronoi point triangle face pt ID
  integer,dimension(:,:),allocatable::vorDel

  ! Voronoi Diagram X,Y,Z coordinates
  real(kind=8),dimension(:),allocatable::vcoordX,vcoordY

  ! Regular structured elements (square cells) with added cells on the edges
  integer,dimension(:,:),allocatable::relmt

  ! Number of grid nodes owned by the processor on a given pet
  integer::sOwnedNode,uOwnedNode

  ! Partition processor ID for unstructured grid
  integer,dimension(:),allocatable::snodeID,selemID,sownedID,snodeLID

  ! Partition processor ID for unstructured grid
  integer,dimension(:),allocatable::unodeID,unodeIDI,uelemID,uownEID,uownedID,unodeLID,delemID

  ! Number of nodes and elements on each partition mesh
  integer::upartN,upartE,spartN,spartE !,upartI

  ! Redistribruted array of arbitrary distrinuted mesh indexes
  integer,dimension(:),allocatable::earthIDs,oceanIDs,spmIDs

  ! Voronoi Cell Declaration Type
  type voronoi
     integer::vertexNb
     integer::border
     integer::btype
     integer::bpoint
     integer,dimension(20)::vertexID
     real(kind=8)::perimeter
     real(kind=8)::area
  end type voronoi
  type(voronoi),dimension(:),allocatable::voronoiCell

  ! Conforming Delaunay Triangulation Declaration Type
  type delaunay
     integer::ngbNb
     integer::sortedHullNb
     integer,dimension(20)::ngbID
     integer,dimension(20)::sortedHull
     real(kind=8)::slope
     real(kind=8),dimension(20)::voronoi_edge
     real(kind=8),dimension(20)::distance
     real(kind=8),dimension(20)::weight
  end type delaunay
  type(delaunay),dimension(:),allocatable::delaunayVertex

  ! Number of seconds per year
  real(kind=8),parameter::secyear=31536000.0

  ! Maximum filling algorithm height
  real(kind=8)::fh

  ! Sediment thickness evolution
  real(kind=8),dimension(:),allocatable::sedchange,tempthick

  ! Infiltration evaporation percentage for water in lakes
  real(kind=8)::infiltration_evaporation

  ! Ice Sheet parameters
  integer::nbix,nbiy,gausskernelSize
  real(kind=8)::ice_dt,ice_dx,ice_xo,ice_yo
  real(kind=8),dimension(:,:),allocatable::iceZ,iceZb,iceH,iceU
  real(kind=8),dimension(:),allocatable::iceX,iceY
  ! Steady time step
  real(kind=8)::ice_Tstep
  ! Deformation constant
  real(kind=8)::ice_deform
  ! Sliding constant
  real(kind=8)::ice_slide
  ! Sliding altitude according to ELA
  real(kind=8)::ice_Zsld
  ! Melting rate below and above ELA
  real(kind=8)::ice_m1,ice_m2

  ! Flexural isostasy parameters
  logical::flexure
  integer::nbfx,nbfy,flex_lay,pressureFields,flexorder

  real(kind=8)::cst1,cst2,cst3,torb
  real(kind=8)::flex_dt,flex_dx,flex_xo,flex_yo,flex_rigid,flex_thick
  real(kind=8)::mean_sediment_density,mean_mantle_density,sea_water_density
  real(kind=8),dimension(:),allocatable::pressTable,poroTable

  real(kind=8),dimension(:,:),allocatable::load,prevload,flexZ,flexDisp,flexSed
  real(kind=8),dimension(:),allocatable::flexX,flexY,sedloader

  ! Ocean parameters
  integer::circON,waveON
  integer::nbox,nboy
  real::storm_dt,ocean_Tstep,latitude,circbase,oceanfric
  real::ocean_dt,ocean_dx,ocean_xo,ocean_yo,courant,oceanfilter
  real,dimension(:,:),allocatable::oceanZ,circZ

  ! Hindcast
  integer::season,forecast
  real(kind=8)::wavebase
  real(kind=8),dimension(:,:),allocatable::ocircU,ocircV,rcircU,rcircV
  real(kind=8),dimension(:,:,:),allocatable::circU,circV

  real,dimension(8)::forecast_param

  ! Number of hindcast scenarios
  type hindcast_param
    ! Percentage of wind/wave subgroup.
    real::oc
    ! Significant wave height (in metres).
    real::hs
    ! Wave period of the energy spectrum
    real::per
    ! Peak wave direction.
    real::dir
    ! Coefficient of directional spreading.
    real::dd
    ! Wind velocity at 10 m elevation (m/s).
    real::wvel
    ! Wind direction.
    real::wdir
  end type hindcast_param

  type hindcast_def
    ! Time start
    real::tstart
    ! Time end
    real::tend
    ! Subgroud parameters
    type(hindcast_param),dimension(:),allocatable::subgroup
  end type hindcast_def
  type(hindcast_def),dimension(:),allocatable::hindcast

  ! Hindcast group number
  integer,dimension(:),allocatable::hcast_gp

  ! Swan model
  character(len=128)::swaninput,swaninfo,swanbot,swanout
  real(kind=8),dimension(:,:,:),allocatable::waveU,waveV
  real(kind=8),dimension(:,:),allocatable::owaveU,owaveV


contains

  ! =====================================================================================
  subroutine term_command(cmds)

    logical(4)::result

    character(len=128)::cmds

    result = .false.

    ! INTEL FORTRAN COMPILER
    ! result = systemqq( cmds )
    ! GNU FORTRAN COMPILER
    call system(cmds)

    return

  end subroutine term_command
  ! =====================================================================================
end module parameters
