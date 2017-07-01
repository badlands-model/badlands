real function SwanIntgratSpc ( p, fmin, fmax, spcsig, theta, wpar, ecs, uloc, vloc, acloc, itype )
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
!   40.87: Marcel Zijlema
!
!   Updates
!
!   40.87, April 2008: New subroutine
!
!   Purpose
!
!   Determine p-th moment of energy density spectrum with respect to a part of frequency space
!   (the integrand may contain power p of relative or absolute frequency, wave number or group velocity)
!
!   Method
!
!   Trapezoidal rule is applied
!
!   Modules used
!
    use ocpcomm4
    use swcomm1
    use swcomm3
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                  :: itype    ! indicate the type of integrand; represent power p of
                                                     ! 1= relative frequency
                                                     ! 2= absolute frequency
                                                     ! 3= wave number
                                                     ! 4= group velocity
    !
    real, intent(in)                     :: p        ! power of the p-th moment
    real, intent(in)                     :: fmax     ! user-specified upper bound of frequency space for integration
    real, intent(in)                     :: fmin     ! user-specified lower bound of frequency space for integration
    real                                 :: uloc     ! ambient u-velocity component at one location
    real                                 :: vloc     ! ambient v-velocity component at one location
    !
    real, dimension(MDC,MSC), intent(in) :: acloc    ! action density at one location
    real, dimension(MDC), intent(in)     :: ecs      ! help array containing (co)sine of spectral directions
    real, dimension(MSC), intent(in)     :: spcsig   ! relative frequency bins
    real, dimension(MDC), intent(in)     :: theta    ! spectral directions
    real, dimension(MSC), intent(in)     :: wpar     ! wave number or group velocity
!
!   Local variables
!
    integer                              :: id       ! loop counter over direction bins
    integer, save                        :: ient = 0 ! number of entries in this subroutine
    integer                              :: is       ! loop counter over frequency bins
    !
    real                                 :: cia      ! first coefficient in trapezoidal rule
    real                                 :: cib      ! second coefficient in trapezoidal rule
    real                                 :: ctail    ! coefficient for high frequency tail in case of (finite) upper bound
    real                                 :: ds       ! width of frequency bin
    real                                 :: dsf      ! remaining width of frequency bin near lower/upper bound
    real                                 :: ead      ! auxiliary factor
    real                                 :: ftail    ! factor for high frequency tail
    real                                 :: ehfr     ! energy at high frequency
    real                                 :: omeg1    ! absolute frequency (upward bin)
    real                                 :: omeg2    ! another absolute frequency (current bin)
    real                                 :: pmom     ! p-th moment of energy density spectrum
    real                                 :: q        ! power of the p-th moment + 1
    real                                 :: us       ! ambient velocity in wave direction
    real                                 :: wloc     ! local wave parameter (wave number or group velocity)
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanIntgratSpc')
    !
    q = p + 1.
    !
    pmom = 0.
    !
    if ( itype == 1 ) then
       !
       do is = 2, MSC
          !
          ds = spcsig(is) - spcsig(is-1)
          !
          if ( spcsig(is-1) >= fmin .and. spcsig(is) <= fmax ) then
             !
             cia = 0.5*spcsig(is-1)**q * ds
             cib = 0.5*spcsig(is  )**q * ds
             !
          elseif ( spcsig(is) > fmax ) then
             !
             dsf = fmax - spcsig(is-1)
             cib = 0.5*fmax**q * dsf**2 / ds
             cia = 0.5*(spcsig(is-1)**q + fmax**q)*dsf - cib
             !
          elseif ( spcsig(is) > fmin ) then
             !
             dsf = spcsig(is) - fmin
             cia = 0.5*fmin**q * dsf**2 / ds
             cib = 0.5*(spcsig(is)**q + fmin**q)*dsf - cia
             !
          else
             !
             cia = 0.
             cib = 0.
             !
          endif
          !
          do id = 1, MDC
             !
             ead  = cia * acloc(id,is-1) + cib * acloc(id,is)
             pmom = pmom + ead * ecs(id) * DDIR
             !
          enddo
          !
          if ( spcsig(is) > fmax ) exit
          !
       enddo
       !
       ! add tail contribution, if appropriate
       !
       if ( fmax > spcsig(MSC) ) then
          !
          if ( MSC > 3 ) then
             !
             if ( PWTAIL(1) <= q ) then
                call msgerr (2, ' power of moment is too large compared to tail power')
                SwanIntgratSpc = -999.
                return
             endif
             !
             ftail = 1. / (PWTAIL(1) - q)
             !
             if ( fmax > 100. ) then
                !
                ctail = 0.
                !
             else
                !
                ctail = fmax**q * (spcsig(MSC)/fmax)**PWTAIL(1)
                !
             endif
             !
             do id = 1, MDC
                !
                ehfr = acloc(id,MSC) * spcsig(MSC) * ecs(id)
                pmom = pmom + ehfr * (spcsig(MSC)**q - ctail) * DDIR * ftail
                !
             enddo
             !
          endif
          !
       endif
       !
    elseif ( itype == 2 ) then
       !
       do id = 1, MDC
          !
          us = uloc*cos(theta(id)+ALCQ) + vloc*sin(theta(id)+ALCQ)
          !
          do is = 2, MSC
             !
             ds = spcsig(is) - spcsig(is-1)
             !
             omeg1 = spcsig(is-1) + wpar(is-1) * us
             omeg2 = spcsig(is  ) + wpar(is  ) * us
             !
             if ( spcsig(is-1) >= fmin .and. spcsig(is) <= fmax ) then
                !
                cia = 0.5*omeg1**p * spcsig(is-1) * ds
                cib = 0.5*omeg2**p * spcsig(is  ) * ds
                !
             elseif ( spcsig(is) > fmax ) then
                !
                dsf   = fmax - spcsig(is-1)
                omeg2 = fmax + (wpar(is)*dsf+wpar(is-1)*(ds-dsf)) * us / ds
                cib   = 0.5*omeg2**p * fmax * dsf**2 / ds
                cia   = 0.5*(omeg1**p * spcsig(is-1) + omeg2**p * fmax)*dsf - cib
                !
             elseif ( spcsig(is) > fmin ) then
                !
                dsf   = spcsig(is) - fmin
                omeg1 = fmin + (wpar(is-1)*dsf+wpar(is)*(ds-dsf)) * us / ds
                cia   = 0.5*omeg1**p * fmin * dsf**2 / ds
                cib   = 0.5*(omeg2**p * spcsig(is) + omeg1**p * fmin)*dsf - cia
                !
             else
                !
                cia = 0.
                cib = 0.
                !
             endif
             !
             ead  = cia * acloc(id,is-1) + cib * acloc(id,is)
             pmom = pmom + ead * ecs(id) * DDIR
             !
             if ( spcsig(is) > fmax ) exit
             !
          enddo
          !
       enddo
       !
       ! add tail contribution, if appropriate
       !
       if ( fmax > spcsig(MSC) ) then
          !
          if ( MSC > 3 ) then
             !
             ftail = 1. / (PWTAIL(1) - 1.)
             !
             if ( fmax > 100. ) then
                !
                ctail = 0.
                !
             else
                !
                ctail = fmax**q * (spcsig(MSC)/fmax)**PWTAIL(1)
                !
             endif
             !
             do id = 1, MDC
                !
                us = uloc*cos(theta(id)+ALCQ) + vloc*sin(theta(id)+ALCQ)
                omeg2 = spcsig(MSC) + wpar(MSC) * us
                !
                ehfr = acloc(id,MSC) * spcsig(MSC) * ecs(id)
                pmom = pmom + ehfr * (omeg2**p * spcsig(MSC) - ctail) * DDIR * ftail
                !
             enddo
             !
          endif
          !
       endif
       !
    elseif ( itype == 3 .or. itype == 4 ) then
       !
       do is = 2, MSC
          !
          ds = spcsig(is) - spcsig(is-1)
          !
          if ( spcsig(is-1) >= fmin .and. spcsig(is) <= fmax ) then
             !
             cia = 0.5*wpar(is-1)**p * spcsig(is-1) * ds
             cib = 0.5*wpar(is  )**p * spcsig(is  ) * ds
             !
          elseif ( spcsig(is) > fmax ) then
             !
             dsf  = fmax - spcsig(is-1)
             wloc = (wpar(is)*dsf+wpar(is-1)*(ds-dsf)) / ds
             cib  = 0.5*wloc**p * fmax * dsf**2 / ds
             cia  = 0.5*(wpar(is-1)**p * spcsig(is-1) + wloc**p * fmax)*dsf - cib
             !
          elseif ( spcsig(is) > fmin ) then
             !
             dsf  = spcsig(is) - fmin
             wloc = (wpar(is-1)*dsf+wpar(is)*(ds-dsf)) / ds
             cia  = 0.5*wloc**p * fmin * dsf**2 / ds
             cib  = 0.5*(wpar(is)**p * spcsig(is) + wloc**p * fmin)*dsf - cia
             !
          else
             !
             cia = 0.
             cib = 0.
             !
          endif
          !
          do id = 1, MDC
             !
             ead  = cia * acloc(id,is-1) + cib * acloc(id,is)
             pmom = pmom + ead * ecs(id) * DDIR
             !
          enddo
          !
          if ( spcsig(is) > fmax ) exit
          !
       enddo
       !
       ! add tail contribution, if appropriate
       !
       if ( fmax > spcsig(MSC) ) then
          !
          if ( MSC > 3 ) then
             !
             ftail = 1. / (PWTAIL(1) - 1.)
             !
             if ( fmax > 100. ) then
                !
                ctail = 0.
                !
             else
                !
                if ( itype == 3 ) then
                   !
                   ctail = (1./GRAV)**p * fmax**(2*p+1) * (spcsig(MSC)/fmax)**PWTAIL(1)
                   !
                elseif ( itype == 4 ) then
                   !
                   ctail = (0.5*GRAV)**p * fmax**(1-p) * (spcsig(MSC)/fmax)**PWTAIL(1)
                   !
                endif
                !
             endif
             !
             do id = 1, MDC
                !
                ehfr = acloc(id,MSC) * spcsig(MSC) * ecs(id)
                pmom = pmom + ehfr * (wpar(MSC)**p * spcsig(MSC) - ctail) * DDIR * ftail
                !
             enddo
             !
          endif
          !
       endif
       !
    endif
    !
    SwanIntgratSpc = pmom
    !
end function SwanIntgratSpc
