!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements flow parameters computation.
module flowcompute

  implicit none

contains

  subroutine discharge(pyStack, pyRcv, pyDischarge, pyDis, pylNodesNb, pygNodesNb)

      integer :: pygNodesNb
      integer :: pylNodesNb
      integer,dimension(pylNodesNb),intent(in) :: pyStack
      integer,dimension(pygNodesNb),intent(in) :: pyRcv
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge

      real(kind=8),dimension(pygNodesNb),intent(out) :: pyDis

      integer :: n, donor, recvr

      pyDis = pyDischarge

      do n = pylNodesNb, 1, -1
        donor = pyStack(n) + 1
        recvr = pyRcv(donor) + 1
        if( donor /= recvr )then
            pyDis(recvr) = pyDis(recvr) + pyDis(donor)
        endif
      enddo

      return

  end subroutine discharge

  subroutine parameters(pyStack, pyRcv, pyDischarge, pyXY, &
      spl_part, pyBid0, pyChi, pyBasinID, pylNodesNb, pygNodesNb)

      integer :: pygNodesNb
      integer :: pylNodesNb
      integer,intent(in) :: pyBid0
      real(kind=8),intent(in) :: spl_part
      integer,dimension(pylNodesNb),intent(in) :: pyStack
      integer,dimension(pygNodesNb),intent(in) :: pyRcv
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
      real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY

      integer,dimension(pygNodesNb),intent(out) :: pyBasinID
      real(kind=8),dimension(pygNodesNb),intent(out) :: pyChi

      integer :: n, donor, recvr, bID
      real(kind=8) :: disch1, disch2, dist

      pyChi = 0.
      pyBasinID = -1
      bID = pyBid0
      do n = 1, pylNodesNb
        donor = pyStack(n) + 1
        recvr = pyRcv(donor) + 1
        if(donor == recvr) bID = bID + 1
        pyBasinID(donor) = bID
        disch1 = pyDischarge(donor)
        disch2 = pyDischarge(recvr)
        if( donor /= recvr .and. disch1 > 0. .and. disch2 > 0.)then
            dist = sqrt( (pyXY(donor,1) - pyXY(recvr,1))**2.0 + &
                (pyXY(donor,2) - pyXY(recvr,2))**2.0 )
            pyChi(donor) = pyChi(recvr) + 0.5*((1./disch2)**spl_part + &
                (1./(disch1))**spl_part) * dist
        endif
      enddo

      return

  end subroutine parameters

  subroutine diffcfl(pyEdges, Cdiff, cfl_dt, pyNodesNb)

      integer :: pyNodesNb
      real(kind=8),intent(in) :: Cdiff
      real(kind=8),dimension(pyNodesNb,20),intent(in) :: pyEdges

      real(kind=8),intent(out) :: cfl_dt

      integer :: p, n
      real(kind=8) :: tmp

      cfl_dt = 1.e6
      do p = 1, pyNodesNb
        n = 1
        do while(pyEdges(p,n) > 0.)
           tmp = 0.05 * pyEdges(p,n)**2. / Cdiff
           cfl_dt = min(tmp,cfl_dt)
           n = n + 1
        enddo
      enddo

      return

  end subroutine diffcfl

  subroutine flowcfl(pyIDs, pyRcv, pyXY, pyElev, pyDischarge, Cero, &
      spl_m, spl_n, cfl_dt, pylNodesNb, pygNodesNb)

      integer :: pygNodesNb
      integer :: pylNodesNb
      real(kind=8),intent(in) :: Cero
      real(kind=8),intent(in) :: spl_m
      real(kind=8),intent(in) :: spl_n
      integer,dimension(pylNodesNb),intent(in) :: pyIDs
      integer,dimension(pygNodesNb),intent(in) :: pyRcv
      real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge

      real(kind=8),intent(out) :: cfl_dt

      integer :: p, d, r
      real(kind=8) :: dz, tmp, dist

      cfl_dt = 1.e6
      do p = 1, pylNodesNb
        d = pyIDs(p) + 1
        r = pyRcv(d) + 1
        dz = pyElev(d) - pyElev(r)
        if( d /= r .and. dz > 0. .and. pyDischarge(d) > 0.)then
            dist = sqrt( (pyXY(d,1)-pyXY(r,1))**2.0 + (pyXY(d,2)-pyXY(r,2))**2.0 )
            tmp = dist / (Cero * pyDischarge(d)**spl_m * (dz/dist)**(spl_n-1.))
            cfl_dt = min(tmp,cfl_dt)
        endif
      enddo

      return

  end subroutine flowcfl

  subroutine sedflux(pyStack, pyRcv, pyXY, pyArea, pyXYmin, pyXYmax, &
      pyMaxH, pyMaxD, pyDischarge, pyFillH, pyElev, pyDiff, Cero, &
      spl_m,spl_n,sea,dt,pyChange,newdt,pylNodesNb,pygNodesNb)

      integer :: pylNodesNb
      integer :: pygNodesNb
      real(kind=8),intent(in) :: dt
      real(kind=8),intent(in) :: sea
      real(kind=8),intent(in) :: spl_n
      real(kind=8),intent(in) :: spl_m
      real(kind=8),intent(in) :: Cero
      real(kind=8),dimension(2),intent(in) :: pyXYmin
      real(kind=8),dimension(2),intent(in) :: pyXYmax
      integer,dimension(pylNodesNb),intent(in) :: pyStack
      integer,dimension(pygNodesNb),intent(in) :: pyRcv
      real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyArea
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyMaxH
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyMaxD
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyFillH
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDiff

      real(kind=8),intent(out) :: newdt
      real(kind=8),dimension(pygNodesNb),intent(out) :: pyChange

      integer :: n, donor, recvr
      real(kind=8) :: maxh, SPL, Qs, dh, mtime, waterH, dist
      real(kind=8),dimension(pygNodesNb) :: sedFluxes

      newdt = dt
      pyChange = -1.e6
      sedFluxes = 0.

      do n = pylNodesNb, 1, -1
        SPL = 0.
        donor = pyStack(n) + 1
        recvr = pyRcv(donor) + 1
        dh = 0.95*(pyElev(donor) - pyElev(recvr))
        if(pyElev(donor) > sea .and. pyElev(recvr) < sea) &
          dh = pyElev(donor) - sea
        if( dh < 0.001 ) dh = 0.
        waterH = pyFillH(donor)-pyElev(donor)

        ! Compute stream power law
        if( recvr /= donor .and. dh > 0.)then
          if(waterH == 0. .and. pyElev(donor) >= sea)then
            dist = sqrt( (pyXY(donor,1)-pyXY(recvr,1))**2.0 + (pyXY(donor,2)-pyXY(recvr,2))**2.0 )
            if(dist > 0.) SPL = -Cero * (pyDischarge(donor))**spl_m * (dh/dist)**spl_n
          endif
        endif

        maxh = pyMaxH(donor)
        if(pyElev(donor) < sea)then
          maxh = sea - pyElev(donor)
        elseif(waterH > 0.)then
          maxh = min( waterH,maxh )
        endif
        maxh = 0.95*maxh
        Qs = 0.

        ! Deposition case
        if( SPL == 0. .and. pyArea(donor) > 0.)then
          if(maxh > 0. .and. pyElev(donor) < sea)then
              if(sedFluxes(donor)*newdt/pyArea(donor) < maxh)then
                SPL = sedFluxes(donor)/pyArea(donor)
                Qs = 0.
              else
                SPL = maxh / newdt
                newdt = min(newdt, maxh/SPL)
                Qs = sedFluxes(donor) - SPL*pyArea(donor)
              endif
          ! Fill depression
          elseif(waterH > 0.0001 .and. donor /= recvr)then
              dh = 0.95*waterH
              if(sedFluxes(donor)*newdt/pyArea(donor) < dh)then
                  SPL = sedFluxes(donor)/pyArea(donor)
                  Qs = 0.
              else
                  SPL = dh / newdt
                  newdt = min(newdt, dh/SPL)
                  Qs = sedFluxes(donor) - SPL*pyArea(donor)
              endif
          ! Base-level (sink)
          elseif(donor == recvr .and. pyArea(donor) > 0.)then
            SPL = sedFluxes(donor) / pyArea(donor)
            Qs = 0.
          else
            Qs = sedFluxes(donor)
          endif
        ! Erosion case
        elseif(SPL < 0.)then
            if(pyElev(donor) > sea .and. pyElev(recvr) < sea) &
                newdt = min( newdt,-0.99*(pyElev(donor)-sea)/SPL)

            if(-SPL * newdt > pyElev(donor) - pyElev(recvr)) &
                newdt = min( newdt,-0.99*(pyElev(donor)-pyElev(recvr))/SPL)

            Qs = -SPL * pyArea(donor) + sedFluxes(donor)
        endif

        ! Update sediment flux in receiver node
        sedFluxes(recvr) = sedFluxes(recvr) + Qs
        pyChange(donor) = SPL + pyDiff(donor)

        ! Update borders
        if(pyXY(donor,1) < pyXYmin(1) .or. pyXY(donor,2) < pyXYmin(2) .or. &
            pyXY(donor,1) > pyXYmax(1) .or. pyXY(donor,2) > pyXYmax(2) ) &
            pyChange(donor) = 0.

        ! Update base levels
        if(donor == recvr .and. pyChange(donor) > 0.)then
            if(waterH == 0.)then
                mtime = pyMaxD(donor) / pyChange(donor)
                newdt = min(newdt, mtime)
            else
                mtime = waterH / pyChange(donor)
                newdt = min(newdt, mtime)
            endif
        endif
      enddo

      return

  end subroutine sedflux

  subroutine sedflux_base(pyStack, pyRcv, pyXY, pyArea, pyXYmin, pyXYmax, &
      pyMaxH, pyDischarge, pyElev, pyDiff, Cero, &
      spl_m,spl_n,sea,dt,pyChange,newdt,pylNodesNb,pygNodesNb)

      integer :: pylNodesNb
      integer :: pygNodesNb
      real(kind=8),intent(in) :: dt
      real(kind=8),intent(in) :: sea
      real(kind=8),intent(in) :: spl_n
      real(kind=8),intent(in) :: spl_m
      real(kind=8),intent(in) :: Cero
      real(kind=8),intent(in) :: pyMaxH
      real(kind=8),dimension(2),intent(in) :: pyXYmin
      real(kind=8),dimension(2),intent(in) :: pyXYmax
      integer,dimension(pylNodesNb),intent(in) :: pyStack
      integer,dimension(pygNodesNb),intent(in) :: pyRcv
      real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyArea
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
      real(kind=8),dimension(pygNodesNb),intent(in) :: pyDiff

      real(kind=8),intent(out) :: newdt
      real(kind=8),dimension(pygNodesNb),intent(out) :: pyChange

      integer :: n, donor, recvr
      real(kind=8) :: maxh, SPL, Qs, dh, mtime, dist
      real(kind=8),dimension(pygNodesNb) :: sedFluxes

      newdt = dt
      pyChange = -1.e6
      sedFluxes = 0.

      do n = pylNodesNb, 1, -1
        SPL = 0.
        donor = pyStack(n) + 1
        recvr = pyRcv(donor) + 1
        dh = 0.95*(pyElev(donor) - pyElev(recvr))
        if(pyElev(donor) > sea .and. pyElev(recvr) < sea) &
          dh = pyElev(donor) - sea
        if(dh < 0.001) dh = 0.

        ! Compute stream power law
        if(recvr /= donor .and. dh > 0.)then
          if(pyElev(donor) >= sea)then
            dist = sqrt( (pyXY(donor,1)-pyXY(recvr,1))**2.0 + (pyXY(donor,2)-pyXY(recvr,2))**2.0 )
            if(dist > 0.) SPL = -Cero * (pyDischarge(donor))**spl_m * (dh/dist)**spl_n
          endif
        endif

        maxh = pyMaxH
        if(pyElev(donor) < sea) maxh = sea - pyElev(donor)
        maxh = 0.99*maxh
        Qs = 0.

        ! Deposition case
        if( SPL == 0. .and. pyArea(donor) > 0.)then
          if(maxh > 0. .and. pyElev(donor) < sea)then
              if(sedFluxes(donor)*newdt/pyArea(donor) < maxh)then
                SPL = sedFluxes(donor)/pyArea(donor)
                Qs = 0.
              else
                SPL = maxh / newdt
                newdt = min(newdt, maxh/SPL)
                Qs = sedFluxes(donor) - SPL*pyArea(donor)
              endif
          ! Base-level (sink)
          elseif(donor == recvr .and. pyArea(donor) > 0.)then
            SPL = sedFluxes(donor) / pyArea(donor)
            Qs = 0.
          else
            Qs = sedFluxes(donor)
          endif
        ! Erosion case
        elseif(SPL < 0.)then
            if(pyElev(donor) > sea .and. pyElev(recvr) < sea) &
                newdt = min( newdt,-0.99*(pyElev(donor)-sea)/SPL)

            if(-SPL * newdt > pyElev(donor) - pyElev(recvr)) &
                newdt = min( newdt,-0.99*(pyElev(donor)-pyElev(recvr))/SPL)

            Qs = -SPL * pyArea(donor) + sedFluxes(donor)
        endif

        ! Update sediment flux in receiver node
        sedFluxes(recvr) = sedFluxes(recvr) + Qs
        pyChange(donor) = SPL + pyDiff(donor)

        ! Update borders
        if(pyXY(donor,1) < pyXYmin(1) .or. pyXY(donor,2) < pyXYmin(2) .or. &
            pyXY(donor,1) > pyXYmax(1) .or. pyXY(donor,2) > pyXYmax(2) ) &
            pyChange(donor) = 0.

        ! Update base levels
        if(donor == recvr .and. pyChange(donor) > 0.)then
            mtime = pyMaxH / pyChange(donor)
            newdt = min(newdt, mtime)
        endif
      enddo

      return

  end subroutine sedflux_base

end module flowcompute
