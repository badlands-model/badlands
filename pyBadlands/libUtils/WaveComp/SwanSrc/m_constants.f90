!------------------------------------------------------------------------------
module m_constants
!------------------------------------------------------------------------------
!
! physical constants
!
real sqrtg   ! square root of grav
real gsq     ! square of grav
real nu      ! kinematic viscosity of water
!
real d_water ! density of water
real d_air   ! density of air
!
real trshdep ! treshold depth (=DEPMIN as given by SWAN)
!
! mathematical constants
!
real pih    ! pi/2
real dera   ! conversion from degrees to radians
real rade   ! conversion from radians to degrees
real expmin ! min argument for exp. function to avoid underflow
real expmax ! max argument for exp. function to avoid overflow
real sqrt2  ! square root of 2 ~ 1.41
!
contains
!
!------------------------------------------------------------------------------
subroutine init_constants
!------------------------------------------------------------------------------
!
use SWCOMM3
!
pih  = 0.5*PI
dera = PI/180.
rade = 180./PI
!
expmin = -20.
expmax =  20.
!
!  physical constants
!
sqrtg   = sqrt(GRAV)
gsq     = GRAV*GRAV
nu      = 1.e-6
d_air   = PWIND(16)
d_water = PWIND(17)
!
trshdep = DEPMIN
!
end subroutine
!
end module
