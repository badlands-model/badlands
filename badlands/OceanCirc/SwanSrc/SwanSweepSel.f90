subroutine SwanSweepSel ( idcmin, idcmax, anybin, iscmin, iscmax, &
                          iddlow, iddtop, idtot , isslow, isstop, &
                          istot , cax   , cay   , rdx   , rdy   , &
                          spcsig)
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
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   computes frequency-dependent counters in directional space
!   and the active bins for considered sweep
!
!   Method
!
!   In order to have stable computation without a CFL restriction
!   certain directional bins lying in the domain of dependence
!   of present vertex are determined. These bins belong to the
!   sweep which is being processed.
!
!   The domain of dependence is enclosed by the two upwave faces
!   of the present vertex in the considered cell. The following
!   set of criterions determine which directional bin belongs to
!   the domain of dependence:
!
!   Cx * rdx(1) + Cy * rdy(1) >= 0 and
!
!   Cx * rdx(2) + Cy * rdy(2) >= 0
!
!   Geometrically, these criterions ensure that the propagation
!   direction towards the present vertex is enclosed between
!   the upwave faces of that vertex in the considered cell.
!
!   The counters of the directional space are frequency dependent.
!   Particularly, the higher frequencies are modified by the ambient
!   current. The lower frequencies (due to the larger wave transport
!   velocity) are less modified by the current.
!
!   Next, it is detemined whether a certain bin lies within a specific
!   sector enclosure of the considered sweep. This is denoted by a
!   logical array anybin.
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(out)                       :: iddlow ! minimum direction bin that is propagated within a sweep
    integer, intent(out)                       :: iddtop ! maximum direction bin that is propagated within a sweep
    integer, intent(out)                       :: idtot  ! maximum number of bins in directional space for considered sweep
    integer, intent(out)                       :: isslow ! minimum frequency that is propagated within a sweep
    integer, intent(out)                       :: isstop ! maximum frequency that is propagated within a sweep
    integer, intent(out)                       :: istot  ! maximum number of bins in frequency space for considered sweep
    !
    integer, dimension(MSC), intent(out)       :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(out)       :: idcmin ! minimum frequency-dependent counter in directional space
    integer, dimension(MDC), intent(out)       :: iscmax ! maximum direction-dependent counter in frequency space
    integer, dimension(MDC), intent(out)       :: iscmin ! minimum direction-dependent counter in frequency space
    !
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cax    ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cay    ! wave transport velocity in y-direction
    real, dimension(2), intent(in)             :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)             :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MSC), intent(in)           :: spcsig ! relative frequency bins
    !
    logical, dimension(MDC,MSC), intent(out)   :: anybin ! true if bin is active in considered sweep
!
!   Local variables
!
    integer       :: id       ! loop counter over direction bins
    integer       :: idclow   ! minimum counter in directional space for given frequency bin
    integer       :: idchgh   ! maximum counter in directional space for given frequency bin
    integer       :: iddum    ! counter in directional space for considered sweep
    integer       :: idsum    ! total active bins in directional space for given frequency bin
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: is       ! loop counter over frequency bins
    integer       :: isclow   ! minimum counter in frequency space for given directional bin
    integer       :: ischgh   ! maximum counter in frequency space for given directional bin
    !
    real          :: caxloc   ! local wave transport velocity in x-direction at present vertex
    real          :: cayloc   ! local wave transport velocity in y-direction at present vertex
    !
    logical       :: lowbin   ! indicates presence of lowest bin in directional space for given frequency bin
    logical       :: hghbin   ! indicates presence of highest bin in directional space for given frequency bin
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanSweepSel')
    !
    ! initialize parameters and arrays
    !
    iddlow =  9999
    iddtop = -9999
    idtot  =     1
    !
    isslow =  9999
    isstop = -9999
    istot  =     1
    !
    idcmin = 0
    idcmax = 0
    anybin = .false.
    !
    iscmin = 1
    iscmax = 1
    !
    ! loop over all frequency bins
    !
    do is = 1, MSC
       !
       ! determine which bin belongs to considered sweep for propagation
       !
       idsum = 0
       !
       if ( is == 1 .or. ICUR /= 0 ) then
          !
          do id = 1, MDC
             !
             caxloc = cax(id,is,1)*rdx(1) + cay(id,is,1)*rdy(1)
             cayloc = cax(id,is,1)*rdx(2) + cay(id,is,1)*rdy(2)
             !
             if ( caxloc >= 0. .and. cayloc >= 0. ) then
                !
                anybin(id,is) = .true.
                idsum         = idsum + 1
                isslow        = min(is,isslow)
                isstop        = max(is,isstop)
                !
             endif
             !
          enddo
          !
       else
          !
          do id = 1, MDC
             !
             ! in case of no current, when first frequency bin is in considered sweep,
             ! other bins for given direction are in this sweep as well
             !
             anybin(id,is) = anybin(id,1)
             !
             if ( anybin(id,1) ) then
                idsum  = idsum + 1
                isstop = max(is,isstop)
             endif
             !
          enddo
          !
       endif
       !
       ! determine sector enclosure in directional space for considered sweep
       !
       idclow = 0
       idchgh = 0
       !
       do id = 1, MDC
          !
          lowbin = .false.
          hghbin = .false.
          !
          if ( anybin(id,is) ) then
             !
             if ( id == 1 ) then
                !
                if ( FULCIR ) then
                   if ( .not.anybin(MDC,is) ) lowbin = .true.
                else
                   lowbin = .true.
                endif
                !
             else
                !
                if ( .not.anybin(id-1,is) ) lowbin = .true.
                !
             endif
             !
             if ( id == MDC ) then
                !
                if ( FULCIR ) then
                   if ( .not.anybin(1,is) ) hghbin = .true.
                else
                  hghbin = .true.
                endif
                !
             else
                !
                if ( .not.anybin(id+1,is) ) hghbin = .true.
                !
             endif
             !
          endif
          !
          if ( lowbin ) idclow = id
          if ( hghbin ) idchgh = id
          !
       enddo
       !
       ! set minimum and maximum counters in directional space for considered sweep
       !
       idcmin(is) = 1
       idcmax(is) = MDC
       !
       if ( idsum == 0 ) then
          !
          idcmin(is) =  9
          idcmax(is) = -9
          !
       elseif ( idsum /= MDC ) then
          !
          if ( idclow > idchgh ) idclow = idclow - MDC
          idcmin(is) = idclow
          idcmax(is) = idchgh
          !
       endif
       !
       if ( idsum /= 0 ) then
          iddlow = min ( iddlow , idcmin(is) )
          iddtop = max ( iddtop , idcmax(is) )
       endif
       !
    enddo
    !
    ! compute maximum number of bins in directional space for considered sweep
    !
    if ( iddlow /= 9999 ) then
       !
       if ( iddtop == -9999 ) then
          call msgerr ( 4, 'inconsistency found in SwanSweepSel: no maximum direction bin ' )
          return
       endif
       !
       idtot = iddtop - iddlow + 1
       !
       if ( ICUR > 0 ) then
          !
          if ( idtot < 3 ) then
             iddtop = iddtop + 1
             if ( idtot == 1 ) iddlow = iddlow - 1
             idtot = 3
          endif
          !
       endif
    else
       !
       if ( iddtop /= -9999 ) then
          call msgerr ( 4, 'inconsistency found in SwanSweepSel: no minimum direction bin ' )
          return
       endif
       !
       idtot = 0
       !
    endif
    !
    if ( idtot > MDC ) then
       iddlow = 1
       iddtop = MDC
       idtot  = MDC
    endif
    !
    ! compute maximum number of bins in frequency space for considered sweep
    !
    if ( isslow /= 9999 ) then
       !
       if ( isstop == -9999 ) then
          call msgerr ( 4, 'inconsistency found in SwanSweepSel: no maximum frequency bin ' )
          return
       endif
       !
       !if ( isslow /= 1 ) then
       !   call msgerr ( 4, 'inconsistency found in SwanSweepSel: isslow <> 1 ' )
       !   return
       !endif
       !
       isslow = 1
       !
       if ( ICUR > 0 ) isstop = max(min(4,MSC),isstop)
       !
       istot = isstop - isslow + 1
       !
    else
       !
       if ( isstop /= -9999 ) then
          call msgerr ( 4, 'inconsistency found in SwanSweepSel: no minimum frequency bin ' )
          return
       endif
       !
       istot = 0
       !
       if ( idtot /= 0 ) then
          call msgerr ( 4, 'inconsistency found in SwanSweepSel: istot = 0 and idtot <> 0 ' )
          return
       endif
       !
    endif
    !
    ! loop over all direction bins
    !
    do iddum = iddlow, iddtop
       id = mod ( iddum - 1 + MDC , MDC ) + 1
       !
       lowbin = .false.
       !
       do is = 1, MSC
          if ( anybin(id,is) ) then
             if ( .not.lowbin ) then
                isclow = is
                lowbin = .true.
             endif
             ischgh = is
          endif
       enddo
       !
       ! set minimum and maximum counters in frequency space for considered sweep
       !
       if ( lowbin ) then
          !
          if ( isclow < isslow .or. ischgh > isstop ) then
             call msgerr ( 4, 'inconsistency found in SwanSweepSel: minimum and maximum counters of frequencies not correct ' )
             return
          endif
          !
          iscmin(id) = isclow
          iscmax(id) = ischgh
          !
       else                     ! no frequency bins fall within considered sweep
          !
          iscmin(id) = 0
          iscmax(id) = 0
          !
       endif
       !
    enddo
    !
end subroutine SwanSweepSel
