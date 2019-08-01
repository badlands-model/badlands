subroutine SwanPropvelX ( cax, cay, ux2, uy2, cgo, ecos, esin )
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
!   41.02: Marcel Zijlema
!
!   Updates
!
!   40.80,     July 2007: New subroutine
!   41.02, February 2009: adaption of velocities in case of diffraction
!
!   Purpose
!
!   computes wave transport velocities of energy in geographical space
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use m_diffr
    use SwanGriddata
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,ICMAX), intent(out) :: cax  ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(out) :: cay  ! wave transport velocity in y-direction
    real, dimension(MSC,ICMAX), intent(in)      :: cgo  ! group velocity
    real, dimension(MDC), intent(in)            :: ecos ! help array containing cosine of spectral directions
    real, dimension(MDC), intent(in)            :: esin ! help array containing sine of spectral directions
    real, dimension(nverts), intent(in)         :: ux2  ! ambient velocity in x-direction at current time level
    real, dimension(nverts), intent(in)         :: uy2  ! ambient velocity in y-direction at current time level
!
!   Local variables
!
    integer       :: ic       ! loop counter over stencil
    integer       :: id       ! loop counter over direction bins
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: is       ! loop counter over frequency bins
    integer       :: ivert    ! vertex index
    !
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanPropvelX')

    do ic = 1, ICMAX
       !
       ivert  = vs(ic)      ! points in computational stencil
       !
       do is = 1, MSC
          do id = 1, MDC
             cax(id,is,ic) = cgo(is,ic) * ecos(id)
             cay(id,is,ic) = cgo(is,ic) * esin(id)
          enddo
       enddo
       !
       ! adapt the celerities in case of diffraction
       !
       if ( IDIFFR /= 0 .and. PDIFFR(3) /= 0. ) then
          do is = 1, MSC
             do id = 1 ,MDC
                cax(id,is,ic) = cax(id,is,ic)*DIFPARAM(ivert)
                cay(id,is,ic) = cay(id,is,ic)*DIFPARAM(ivert)
             enddo
          enddo
       endif
       !
       ! ambient currents added
       !
       if ( ICUR /= 0 )  then
          do is = 1, MSC
             do id = 1, MDC
                cax(id,is,ic) = cax(id,is,ic) + ux2(ivert)
                cay(id,is,ic) = cay(id,is,ic) + uy2(ivert)
             enddo
          enddo
       endif
       !
    enddo
    !
end subroutine SwanPropvelX
