!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the pyReef carbonate platform modelling application.     !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Main wave modelling data class.

module classdata

  ! SP grid model values
  integer::sp_n,sp_m
  integer::xlen,ylen,xclen,yclen,stratal_x,stratal_y
  real::stratal_xo,stratal_yo,stratal_dx,stratal_xm,stratal_ym
  real::sea_level
  real,dimension(:,:),allocatable::sp_topo

  ! Forecast parameters
  real::wave_base
  real,dimension(6)::forecast_param

  ! Number of hindcast scenarios
  type hindcast_param
    ! Wave height.
    real::wh
    ! Wave period.
    real::wp
    ! Wave/wind direction.
    real::wdir
    ! Wind velocity at 10 m elevation (m/s).
    real::wvel
    ! Wind direction.
    real::wddir
  end type hindcast_param
  type(hindcast_param)::hindcast

end module classdata
! ============================================================================
module mpidata

#include "mpif.h"

  ! Error code
  integer::ierr
  ! Processor ID
  integer::iam
  ! Number of processors
  integer::nprocs

end module mpidata
! ============================================================================
module miscdata

  ! Input / output files
  character(len=128)::finput
  character(len=128)::xyzfile
  character(len=128)::swaninput
  character(len=128)::swaninfo
  character(len=128)::swanbot
  character(len=128)::swanout
  character(len=128)::h5data

  ! MPI integer type communicator
  integer::int_type
  ! MPI double type communicator
  integer::dbl_type
  ! MPI double type communicator
  integer::real_type
  ! MPI logical type communicator
  integer::lgc_type
  ! MPI max type communicator
  integer::max_type
  ! MPI min type communicator
  integer::min_type
  ! MPI sum type communicator
  integer::sum_type
  ! SPModel communicator
  integer::ocean_comm_world

contains

  ! ============================================================================
  subroutine append_str(stg1,stg2)

    integer::l1,l2
    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2

    return

  end subroutine append_str
  ! ============================================================================
  subroutine noblnk(string)

    integer::i,j,lg
    character(len=128)::string

    lg=len(string)
    do
      if(lg<=0.or.string(lg:lg)/=' ') exit
      lg=lg-1
    enddo
    if(lg>0)then
      ! find first non-blank character
      i=1
      do
        if(i>lg.or.(string(i:i)/=' '.and.string /= ' ')) exit
        i=i+1
      enddo
      ! determine end of continuous (non-blank) string
      j=i
      do
        if(j>lg)then
          exit
        elseif(string(j:j)==' ')then
          exit
        elseif(string=='  ')then
          exit
        elseif(j==128)then
          exit
        endif
        j=j+1
      enddo
      ! j points to first blank position or past end of string; adjust to last
      ! non-blank position in string
      j=min(j-1,lg)
      string=string(i:j)
      if(j<len(string)) string(j+1:len(string))=' '
    else
       ! there were only blanks in string
       string=' '
    endif

    return

  end subroutine noblnk
  ! ============================================================================

end module miscdata
! ============================================================================
