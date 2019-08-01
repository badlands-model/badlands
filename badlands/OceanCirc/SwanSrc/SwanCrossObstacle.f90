logical function SwanCrossObstacle ( xv, yv, xobs, yobs )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2011  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.80: Marcel Zijlema
!
!   Updates
!
!   40.80, March 2008: New subroutine
!
!   Purpose
!
!   Checks whether an obstacle line cross a face in computational grid
!
!   Method
!
!   See Technical documentation
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Argument variables
!
    real, dimension(2), intent(in) :: xobs     ! x-coordinate of obstacle point
    real, dimension(2), intent(in) :: xv       ! x-coordinate of vertex of face
    real, dimension(2), intent(in) :: yobs     ! y-coordinate of obstacle point
    real, dimension(2), intent(in) :: yv       ! y-coordinate of vertex of face
!
!   Local variables
!
    integer, save                  :: ient = 0 ! number of entries in this subroutine
    !
    real                           :: a        ! dummy variable
    real                           :: b        ! dummy variable
    real                           :: c        ! dummy variable
    real                           :: d        ! dummy variable
    real                           :: det      ! determinant
    real                           :: e        ! dummy variable
    real                           :: f        ! dummy variable
    real                           :: p        ! dummy variable
    real                           :: q        ! dummy variable
    !
    logical                        :: EQREAL   ! indicate whether two reals are equal or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanCrossObstacle')
    !
    ! initially, we assume there is crossing
    !
    SwanCrossObstacle = .true.
    !
    a = xv(1) - xobs(1)
    b = yv(1) - yobs(1)
    !
    c = xv(2) - xv(1)
    d = yv(2) - yv(1)
    !
    e = xobs(2) - xobs(1)
    f = yobs(2) - yobs(1)
    !
    ! compute determinant
    !
    det = e*d - f*c
    !
    if ( .not.EQREAL(det,0.) ) then
       !
       p = (a*d - b*c)/det
       q = (a*f - b*e)/det
       !
       if ( p<0. .or. p>1. .or. q<0. .or. q>1. ) SwanCrossObstacle = .false.
       !
    else
       !
       SwanCrossObstacle = .false.
       !
    endif
    !
end function SwanCrossObstacle
