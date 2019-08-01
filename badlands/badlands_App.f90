! =====================================================================================
! =====================================================================================
! BADLANDS (BAsin and LAndscape Dynamics)
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
!       Filename:  badlands_App.f90
!
!    Description:  Top level BADLANDS Application Driver
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
program BADLANDS_Application

  use restart
  use geomesh
  use coupling
  use parallel
  use geomorpho
  use stratmorph
  use parameters
  use strata_evol

  implicit none

  integer::opt
  real(kind=8)::t1,t2

  ! start up MPI
  call mpi_init(rc)
  call mpi_comm_size(badlands_world,npets,rc)
  call mpi_comm_rank(badlands_world,pet_id,rc)

  ! Get experiment file name
  xmlfile=''
  opt=iargc()
  if(opt<1)then
    print*,'Wrong command line, use: mpirun -n X ./badlands <xml-file-name>'
    call mpi_finalize(rc)
  elseif(opt==1)then
    call getarg(1,xmlfile)
  endif

  ! Define simulation meshes
  t1=mpi_wtime()
  call GeoMesher
  call StrataGen
  t2=mpi_wtime()
  if(pet_id==0)print*,'-------------------------'
  if(pet_id==0)print*,'BADLANDS Grid Components Initialized (s) ',t2-t1
  if(pet_id==0)print*,'-------------------------'

  ! Interpolation of regular grid component to unstructured one
  t1=mpi_wtime()
  if(restartFlag)then
    if(.not.disp3d)then
      call getSPM_hdf5topography
    else
      call getRestart_topography
    endif
  else
    call bilinearTopo
  endif
  call bilinearGrid
  t2=mpi_wtime()
  if(pet_id==0)print*,'-------------------------'
  if(pet_id==0)print*,'BADLANDS Bilinear Interpolation Initialized (s) ',t2-t1
  if(pet_id==0)print*,'-------------------------'

  do while(simulation_time<time_end)

    ! Perform ice sheet flow modelling
    if(ice_dx>0.) call getIceModel

    ! Perform ocean and wave velocity modelling
    if(ocean_dx>0.) call getOceanModel

    ! Perform ice sheet flow modelling
    if(flexure) call getFlexModel

    ! Get the current topographic state, performs Geodynamic Evolution calculation and
    ! prepare export displacement field arrays for the SPM model
    if(simulation_time>time_start) call getEarthData
    ! If the displacement field has a horizontal component change
    if(disp3d)then
      t1=mpi_wtime()
      call mvSpmGrid
      t2=mpi_wtime()
      if(pet_id==0)print*,'-------------------------'
      if(pet_id==0)print*,'BADLANDS Grid Horizontal Displacement (s) ',t2-t1
      if(pet_id==0)print*,'-------------------------'
    endif

    ! Based on precipitation rate and displacement field compute landscape and geomorphological
    ! evolution using the SPM model
    t1=mpi_wtime()
    if(totgrn>0)then
      call stratgeomorph
    else
      call geomorphology
    endif
    t2=mpi_wtime()
    if(pet_id==0)print*,'-------------------------'
    if(pet_id==0)print*,'BADLANDS Surface Process Model (s) ',t2-t1
    if(pet_id==0)print*,'-------------------------'
  enddo

  call mpi_finalize(rc)

end program BADLANDS_Application
