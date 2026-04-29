!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Entry call to compute the precipitation using the linear orographic precipitation
! model of Smith and Barstad (2004).
subroutine compute( pyElev, pyDx, pywindX, pywindY, pyminRain, pymaxRain, pybackRain, pyNm, &
  pyCw, pyHw, pytauC, pytauF, pyRain, pynx, pyny)

  use classoro
  implicit none

  integer :: pynx
  integer :: pyny
  real,intent(in) :: pyDx
  real,intent(in) :: pywindX
  real,intent(in) :: pywindY
  real,intent(in) :: pyminRain
  real,intent(in) :: pymaxRain
  real,intent(in) :: pybackRain
  real,intent(in) :: pyNm
  real,intent(in) :: pyCw
  real,intent(in) :: pyHw
  real,intent(in) :: pytauC
  real,intent(in) :: pytauF
  real,intent(in) :: pyElev(pynx,pyny)
  real,intent(out) :: pyRain(pynx,pyny)

  integer :: i, j, p, q
  integer :: i1, j1, i2, j2

  ! Conversion time -- 200 to 2000 s
  tauc = pytauC
  ! Fallout time -- 200 to 2000 s
  tauf = pytauF
  ! Background precipitation rate -- 0 to 5 m/a
  prbgd = pybackRain
  ! Rain min/max range
  minRain = pyminRain
  maxRain = pymaxRain
  ! Horizontal wind components (m/s) -- 1 to 100 m/s
  u = pywindX
  v = pywindY
  ! Moist stability frequency -- 0 to 0.01 s 1
  nm = pyNm
  ! Uplift sensitivity factor -- 0.001 to 0.02 kg/m3
  cw = pyCw
  ! Water vapor scale height -- 1 to 5 km in m
  hw = pyHw
  ! Grid size
  dx = pyDx

  ! Define computational grid
  nx = int(pynx/2.) + pynx
  ny = int(pyny/2.) + pyny
  i1 = int(pynx/4.)
  j1 = int(pyny/4.)
  i2 = int(pynx/4.) + pynx - 1
  j2 = int(pyny/4.) + pyny - 1

  ! Read elevation file
  if(allocated(elev)) deallocate(elev)
  allocate(elev(nx,ny))

  elev(:,:) = 0.
  q = 1
  do j = j1, j2
    p = 1
    do i = i1, i2
      elev(i,j) = pyElev(p,q)
      p = p + 1
    enddo
    q = q + 1
  enddo

  ! Extrapolate on borders
  call extrapolate_border(i1, j1, i2, j2)

  ! Compute orographic rain
  call get_rain

  ! Define precipitation using user range
  if(allocated(rainval)) deallocate(rainval)
  allocate(rainval(pynx,pyny))

  call clip_rain(i1, j1, i2, j2, pynx, pyny)

  ! Allocate orographic rain value
  do j = 1, pyny
    do i = 1, pynx
      pyRain(i,j) = rainval(i,j)
    enddo
  enddo

  return

end subroutine compute
