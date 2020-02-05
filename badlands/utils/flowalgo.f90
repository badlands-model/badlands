!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Loop over node coordinates and find if they belong to local partition.
subroutine overlap(pyX,pyY,pyXst,pyYst,pyXed,pyYed,pyPart,pyNodes)
  implicit none

  integer :: q
  integer :: pyNodes
  real(kind=8),intent(in) :: pyXst
  real(kind=8),intent(in) :: pyXed
  real(kind=8),intent(in) :: pyYst
  real(kind=8),intent(in) :: pyYed
  real(kind=8),intent(in) :: pyX(pyNodes)
  real(kind=8),intent(in) :: pyY(pyNodes)
  integer,intent(out) :: pyPart(pyNodes)

  pyPart = -1
  do q = 1,pyNodes
      if( pyX(q) >= pyXst .and. pyX(q) <= pyXed )then
          if( pyY(q) >= pyYst .and. pyY(q) <= pyYed )then
              pyPart(q) = q-1
          endif
      endif
  enddo

  return

end subroutine overlap

! This module builds the stack used for in Braun & Willett (2013) flow network
! algorithm
subroutine build(pyBase,pyRcv,pyDelta,pyDonors,pyStackOrder,pyBaseNb,pyDeltaNb,pyNodesNb)

  use classfv
  implicit none

  integer :: pyBaseNb
  integer :: pyDeltaNb
  integer :: pyNodesNb
  integer,dimension(pyBaseNb),intent(in) :: pyBase
  integer,dimension(pyNodesNb),intent(in) :: pyRcv
  integer,dimension(pyDeltaNb),intent(in) :: pyDelta

  integer,dimension(pyNodesNb),intent(out) :: pyDonors
  integer,dimension(pyNodesNb),intent(out) :: pyStackOrder

  integer :: p,j,k,success
  integer,dimension(pyNodesNb) :: intArray

  j = 0

  if(allocated(stackOrder)) deallocate(stackOrder)
  if(allocated(allocs)) deallocate(allocs)
  if(allocated(Donors)) deallocate(Donors)
  if(allocated(Delta)) deallocate(Delta)
  allocate(stackOrder(pyNodesNb))
  allocate(allocs(pyNodesNb))
  allocate(Donors(pyNodesNb))
  allocate(Delta(pyDeltaNb))

  intArray = 0

  stackOrder = 0
  Delta = pyDelta+1

  do k = 1, pyNodesNb
      Donors(Delta(pyRcv(k)+1) + intArray(pyRcv(k)+1)) = k
      intArray(pyRcv(k)+1) = intArray(pyRcv(k)+1)+1
  enddo

  allocs = -1
  do p = 1, pyBaseNb
      k = pyBase(p)+1
      j = j+1
      stackOrder(j) = k
      allocs(k) = p
      success = addtostack(p,k,j)
  enddo

  pyDonors = Donors-1
  pyStackOrder = stackOrder-1

  return

end subroutine build

! This module implements flow parameters computation.
subroutine eroparams(typefct, m, n, mt, nt, kt, kw, b, bsfct)

  use classfv
  implicit none

  integer :: typefct
  integer :: bsfct
  real(kind=8),intent(in) :: m
  real(kind=8),intent(in) :: n
  real(kind=8),intent(in) :: mt
  real(kind=8),intent(in) :: nt
  real(kind=8),intent(in) :: kt
  real(kind=8),intent(in) :: kw
  real(kind=8),intent(in) :: b

  incisiontype = typefct
  bedslptype = bsfct
  spl_m = m
  spl_n = n
  sed_mt = mt
  sed_nt = nt
  sed_kt = kt
  width_kw = kw
  width_b = b

  return

end subroutine eroparams

subroutine discharge(sea, pyStack, pyRcv, pyElev, pyDischarge, pyDis, pyLay, pylNodesNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pygNodesNb
  integer :: pylNodesNb
  integer,dimension(pylNodesNb),intent(in) :: pyStack
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
  real(kind=8),intent(in) :: sea
  real(kind=8),dimension(pygNodesNb),intent(out) :: pyDis
  real(kind=8),dimension(pygNodesNb),intent(out) :: pyLay

  integer :: n, donor, recvr

  pyDis = pyDischarge
  pyLay = 0.

  do n = pylNodesNb, 1, -1
    donor = pyStack(n) + 1
    recvr = pyRcv(donor) + 1
    if( donor /= recvr )then
  ! does not sum discharge below sea level - only one donor per receiver
        if (pyElev(recvr)>=sea)then
            pyDis(recvr) = pyDis(recvr) + pyDis(donor)
        else
            if (pyDis(donor)>=pyDis(recvr)) then
              pyDis(recvr)=pyDis(donor)
            endif
        endif
    endif
    pyLay(donor) = pyElev(donor)-pyElev(recvr)
  enddo

  return

end subroutine discharge

subroutine parameters(pyStack, pyRcv, pyDischarge, pyXY, pyBid0, pyChi, pyBasinID, pylNodesNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pygNodesNb
  integer :: pylNodesNb
  integer,intent(in) :: pyBid0
  integer,dimension(pylNodesNb),intent(in) :: pyStack
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
  real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY

  integer,dimension(pygNodesNb),intent(out) :: pyBasinID
  real(kind=8),dimension(pygNodesNb),intent(out) :: pyChi

  integer :: n, donor, recvr, bID
  real(kind=8) :: disch1, disch2, dist, slp2

  if(spl_n > 0)then
    slp2 = spl_m / spl_n
  else
    slp2 = 1.
  endif
  pyChi = 0.
  pyBasinID = -1
  bID = pyBid0
  do n = 1, pylNodesNb
    donor = pyStack(n) + 1
    recvr = pyRcv(donor) + 1
    if(recvr == 0) donor = recvr
    if(donor == recvr) bID = bID + 1
    if(donor>0)then
      pyBasinID(donor) = bID
      disch1 = pyDischarge(donor)
      disch2 = pyDischarge(recvr)
      if( donor /= recvr .and. disch1 > 0. .and. disch2 > 0.)then
          dist = sqrt( (pyXY(donor,1) - pyXY(recvr,1))**2.0 + &
              (pyXY(donor,2) - pyXY(recvr,2))**2.0 )
          if(spl_n > 0)then
            pyChi(donor) = pyChi(recvr) + 0.5*((1./disch2)**slp2 + &
              (1./(disch1))**slp2) * dist
          else
            pyChi(donor) = 0.
          endif
      endif
    endif
  enddo

  return

end subroutine parameters

subroutine basinparameters(pyStack, pyRcv, pyElev, pyWatH, pyArea, &
pyBasinID, pyVolume, pylNodesNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pygNodesNb
  integer :: pylNodesNb
  integer,dimension(pylNodesNb),intent(in) :: pyStack
  integer,dimension(pygNodesNb),intent(in) :: pyRcv

  real(kind=8),dimension(pygNodesNb),intent(in) :: pyWatH
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyArea

  integer,dimension(pygNodesNb),intent(out) :: pyBasinID
  real(kind=8),dimension(pygNodesNb),intent(out) :: pyVolume

  integer :: n, donor, recvr, pitID

  pyBasinID = -1
  pyVolume = -1
  pitID = -1
  do n = 1, pylNodesNb
    donor = pyStack(n) + 1
    recvr = pyRcv(donor) + 1
    if(donor == recvr)then
      pitID = donor - 1
      pyBasinID(donor) = pitID
    endif
    if(pyWatH(donor) > pyElev(donor) .and. pitID > -1)then
      pyBasinID(donor) = pitID
      pyVolume(pitID+1) = pyVolume(pitID+1) + (pyWatH(donor)-pyElev(donor))*pyArea(donor)
    endif
  enddo

  return

end subroutine basinparameters

subroutine basindrainage(orderPits, pitID, pyRcv, pIDs, fillH, &
sea, pyDrain, pitNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pitNb
  integer :: pygNodesNb
  real(kind=8),intent(in) :: sea
  integer,dimension(pitNb),intent(in) :: orderPits
  integer,dimension(pygNodesNb),intent(in) :: pitID
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  integer,dimension(pitNb),intent(in) :: pIDs

  real(kind=8),dimension(pygNodesNb),intent(in) :: fillH

  integer,dimension(pygNodesNb),intent(out) :: pyDrain

  integer,dimension(pitNb) :: chainDrain
  integer :: n, donor, recvr, nID, count, p, newDrain
  logical :: newpit,exist

  pyDrain = -1
  chainDrain = -1
  count = 1

  do n = 1, pitNb
    nID = pIDs(orderPits(n)+1) + 1
    donor = nID
    exist = .False.
    if(pyDrain(nID)>-1) exist = .True.
    if(fillH(nID)<sea)then
      exist = .True.
      pyDrain(nID) = nID-1
    endif
    do while(.not. exist)
      recvr = pyRcv(donor) + 1
      ! If this is an internal drained basin or an edge node
      if(recvr == donor)then
        pyDrain(nID) = recvr - 1
        count = 1
        chainDrain = -1
        exist = .True.
      elseif(fillH(recvr)<sea)then
        pyDrain(nID) = recvr - 1
        count = 1
        chainDrain = -1
        exist = .True.
      elseif(pitID(recvr) == -1 .or. pitID(recvr) == nID-1)then
        donor = recvr
      elseif(fillH(pitID(recvr))<sea)then
        donor = recvr
      else
        p = 1
        newpit = .True.
        newDrain = pitID(recvr)
        do while(chainDrain(p) >= 0)
           if(chainDrain(p)==newDrain) newpit = .False.
           p = p + 1
        enddo
        if(newpit)then
          chainDrain(p) = newDrain
          pyDrain(nID) = newDrain
          donor = newDrain + 1
          nID = newDrain + 1
        else
          donor = recvr
        endif
      endif
    enddo
  enddo

  return

end subroutine basindrainage

subroutine basindrainageall(orderPits, pitID, pyRcv, pIDs, &
pyDrain, pitNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pitNb
  integer :: pygNodesNb
  integer,dimension(pitNb),intent(in) :: orderPits
  integer,dimension(pygNodesNb),intent(in) :: pitID
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  integer,dimension(pitNb),intent(in) :: pIDs

  integer,dimension(pygNodesNb),intent(out) :: pyDrain

  integer,dimension(pitNb) :: chainDrain
  integer :: n, donor, recvr, nID, count, p, newDrain
  logical :: newpit,exist

  pyDrain = -1
  chainDrain = -1
  count = 1

  do n = 1, pitNb
    nID = pIDs(orderPits(n)+1) + 1
    donor = nID
    exist = .False.
    if(pyDrain(nID)>-1) exist = .True.
    do while(.not. exist)
      recvr = pyRcv(donor) + 1
      ! If this is an internal drained basin or an edge node
      if(recvr == donor)then
        pyDrain(nID) = recvr - 1
        count = 1
        chainDrain = -1
        exist = .True.
      elseif(pitID(recvr) == -1 .or. pitID(recvr) == nID-1)then
        donor = recvr
      else
        p = 1
        newpit = .True.
        newDrain = pitID(recvr)
        do while(chainDrain(p) >= 0)
           if(chainDrain(p)==newDrain) newpit = .False.
           p = p + 1
        enddo
        if(newpit)then
          chainDrain(p) = newDrain
          pyDrain(nID) = newDrain
          donor = newDrain + 1
          nID = newDrain + 1
        else
          donor = recvr
        endif
      endif
    enddo
  enddo

  return

end subroutine basindrainageall

subroutine flowcfl(pyIDs, pyRcv, pyXY, pyElev, pyDischarge, Cero, &
cfl_dt, pylNodesNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pygNodesNb
  integer :: pylNodesNb
  integer,dimension(pylNodesNb),intent(in) :: pyIDs
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  real(kind=8),dimension(pygNodesNb),intent(in) :: Cero
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
        tmp = dist / (Cero(d) * pyDischarge(d)**spl_m * (dz/dist)**(spl_n-1.))
        cfl_dt = min(tmp,cfl_dt)
    endif
  enddo

  return

end subroutine flowcfl

subroutine diffmarine(pyZ, pyBord, pyDepoH, pyNgbs, pyEdge, pyDist, pyCoeff, pyGIDs, &
                  slvl, pymaxth, tstep, pyDiff, mindt, pylocalNb, pyglobalNb)

  use classfv
  implicit none

  integer :: pyglobalNb
  integer :: pylocalNb
  integer,dimension(pylocalNb),intent(in) :: pyGIDs
  integer,dimension(pyglobalNb),intent(in) :: pyBord
  integer,dimension(pyglobalNb,20),intent(in) :: pyNgbs

  real(kind=8),intent(in) :: slvl
  real(kind=8),intent(in) :: pymaxth
  real(kind=8),intent(in) :: tstep
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyCoeff
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyDepoH
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist

  real(kind=8),intent(out) :: mindt
  real(kind=8),dimension(pyglobalNb),intent(out) :: pyDiff

  integer :: k, gid, ngbid, p
  real(kind=8) :: flx

  pyDiff = 0.
  mindt = tstep
  do k = 1, pylocalNb
    gid = pyGIDs(k)+1
    if(pyBord(gid)>0 .and. pyZ(gid)<slvl)then
      loop: do p =1,20
        if(pyNgbs(gid,p)<0) exit loop
        ngbid = pyNgbs(gid,p)+1
        if(pyBord(ngbid)>0.)then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          elseif(pyDepoH(ngbid)>pymaxth .and. pyZ(gid)<pyZ(ngbid) .and. pyZ(ngbid)<slvl)then
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          endif
        elseif(pyBord(ngbid)<1)then
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          endif
        endif
      enddo loop
      ! In case we diffuse more sediment than what is available during the
      ! considered time step, flag it
      if(pyDiff(gid)<0. .and. pyDiff(gid)*tstep<-pyDepoH(gid))then
        mindt = min(-pyDepoH(gid)/pyDiff(gid),mindt)
      endif
    endif
  enddo

  return

end subroutine diffmarine

subroutine difffailure(pyZ, pyBord, pyDepoH, pyNgbs, pyEdge, pyDist, pyCoeff, pyGIDs, &
                  pymaxth, tstep, pyDiff, mindt, pylocalNb, pyglobalNb)

  use classfv
  implicit none

  integer :: pyglobalNb
  integer :: pylocalNb
  integer,dimension(pylocalNb),intent(in) :: pyGIDs
  integer,dimension(pyglobalNb),intent(in) :: pyBord
  integer,dimension(pyglobalNb,20),intent(in) :: pyNgbs

  real(kind=8),intent(in) :: pymaxth
  real(kind=8),intent(in) :: tstep
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyCoeff
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyDepoH
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist

  real(kind=8),intent(out) :: mindt
  real(kind=8),dimension(pyglobalNb),intent(out) :: pyDiff

  integer :: k, gid, ngbid, p
  real(kind=8) :: flx

  pyDiff = 0.
  mindt = tstep
  do k = 1, pylocalNb
    gid = pyGIDs(k)+1
    if(pyBord(gid)>0)then
      loop: do p =1,20
        if(pyNgbs(gid,p)<0) exit loop
        ngbid = pyNgbs(gid,p)+1
        if(pyBord(ngbid)>0.)then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          elseif(pyDepoH(ngbid)>pymaxth .and. pyZ(gid)<pyZ(ngbid))then
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          endif
        elseif(pyBord(ngbid)<1)then
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
            pyDiff(gid) = pyDiff(gid) + pyCoeff(gid)*flx
          endif
        endif
      enddo loop
      ! In case we diffuse more sediment than what is available during the
      ! considered time step, flag it
      if(pyDiff(gid)<0. .and. pyDiff(gid)*tstep<-pyDepoH(gid))then
        mindt = min(-pyDepoH(gid)/pyDiff(gid),mindt)
      endif
    endif
  enddo

  return

end subroutine difffailure

subroutine diffsedmarine(pyZ, pyBord, pyDepo, pyDepoH, slvl, pymaxth, pyCoeff, pyNgbs, pyEdge, &
                    pyDist, pyGIDs, pyDiff, sumDiff, pylocalNb, pyglobalNb, pyRockNb)

  use classfv
  implicit none

  integer :: pyglobalNb
  integer :: pylocalNb
  integer :: pyRockNb
  integer,dimension(pylocalNb),intent(in) :: pyGIDs
  integer,dimension(pyglobalNb),intent(in) :: pyBord
  integer,dimension(pylocalNb,20),intent(in) :: pyNgbs

  real(kind=8),intent(in) :: slvl
  real(kind=8),intent(in) :: pymaxth
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyCoeff
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyDepoH
  real(kind=8),dimension(pyglobalNb,pyRockNb),intent(in) :: pyDepo

  real(kind=8),dimension(pyglobalNb,pyRockNb),intent(out) :: pyDiff
  real(kind=8),dimension(pyglobalNb),intent(out) :: sumDiff

  integer :: k, gid, ngbid, p, r
  real(kind=8) :: flx, frac, sfrac, tfrac, sed(pyRockNb), tsed

  pyDiff = 0.
  sumDiff = 0.

  do k = 1, pylocalNb
    gid = pyGIDs(k)+1
    if(pyBord(gid)>0 .and. pyZ(gid)<slvl)then
      loop: do p =1,20
        if(pyNgbs(gid,p)<0) exit loop
        ngbid = pyNgbs(gid,p)+1
        if(pyBord(ngbid)>0.)then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            sfrac = 0.
            sed = 0.
            tsed = 0.
            do r = 1, pyRockNb
              frac = pyDepo(gid,r)/pyDepoH(gid)
              sfrac = sfrac + frac
              sed(r) = pyCoeff(gid)*frac*flx
              tsed = tsed + sed(r)
            enddo
            if(sfrac>0.)then
              tfrac = 1./sfrac
              pyDiff(gid,:) = pyDiff(gid,:) + tfrac*sed(:)
              sumDiff(gid) = sumDiff(gid) + tfrac*tsed
            endif
          elseif(pyDepoH(ngbid)>pymaxth .and. pyZ(gid)<pyZ(ngbid) .and. pyZ(ngbid)<slvl)then
            sfrac = 0.
            sed = 0.
            tsed = 0.
            do r = 1, pyRockNb
              frac = pyDepo(ngbid,r)/pyDepoH(ngbid)
              sfrac = sfrac + frac
              sed(r) = pyCoeff(gid)*frac*flx
              tsed = tsed + sed(r)
            enddo
            if(sfrac>0.)then
              tfrac = 1./sfrac
              pyDiff(gid,:) = pyDiff(gid,:) + tfrac*sed(:)
              sumDiff(gid) = sumDiff(gid) + tfrac*tsed
            endif
          endif
        elseif(pyBord(ngbid)<1)then
          if(pyDepoH(gid)>pymaxth .and. pyZ(gid)>pyZ(ngbid))then
            sfrac = 0.
            sed = 0.
            tsed = 0.
            flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
            do r = 1, pyRockNb
              frac = pyDepo(gid,r)/pyDepoH(gid)
              sfrac = sfrac + frac
              sed(r) = pyCoeff(gid)*frac*flx
              tsed = tsed + sed(r)
            enddo
            if(sfrac>0.)then
              tfrac = 1./sfrac
              pyDiff(gid,:) = pyDiff(gid,:) + tfrac*sed(:)
              sumDiff(gid) = sumDiff(gid) + tfrac*tsed
            endif
          endif
        endif
      enddo loop
    endif
  enddo

  return

end subroutine diffsedmarine

subroutine diffsedhillslope(pyZ, pyBord, difflay, maxlayh, pyCoeff, pyNgbs, pyEdge, pyDist, pyGIDs,  &
                     sumDiff, ero, depo, pylocalNb, pyglobalNb, pyRockNb)

  use classfv
  implicit none

  integer :: pyglobalNb
  integer :: pylocalNb
  integer :: pyRockNb
  integer,dimension(pylocalNb),intent(in) :: pyGIDs
  integer,dimension(pyglobalNb),intent(in) :: pyBord
  integer,dimension(pylocalNb,20),intent(in) :: pyNgbs

  real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
  real(kind=8),dimension(pyglobalNb),intent(in) :: pyCoeff
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
  real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist
  real(kind=8),dimension(pyglobalNb),intent(in) :: maxlayh
  real(kind=8),dimension(pyglobalNb,pyRockNb),intent(in) :: difflay

  real(kind=8),dimension(pyglobalNb),intent(out) :: sumDiff
  real(kind=8),dimension(pyglobalNb,pyRockNb),intent(out) :: depo
  real(kind=8),dimension(pyglobalNb,pyRockNb),intent(out) :: ero

  integer :: k, gid, ngbid, p, r
  real(kind=8) :: flx, frac, sfrac, tfrac, sed(pyRockNb), tsed

  ero = 0.
  depo = 0.
  sumDiff = 0.

  do k = 1, pylocalNb
  gid = pyGIDs(k)+1
  if(pyBord(gid)>0)then
    loop: do p =1,20
      if(pyNgbs(gid,p)<0) exit loop
      ngbid = pyNgbs(gid,p)+1
      if(pyBord(ngbid)>0.)then
        if(pyZ(gid)>pyZ(ngbid))then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          sfrac = 0.
          sed = 0.
          tsed = 0.
          do r = 1, pyRockNb
            frac = difflay(gid,r)/maxlayh(gid)
            sfrac = sfrac + frac
            sed(r) = pyCoeff(gid)*frac*flx
            tsed = tsed + sed(r)
          enddo
          if(sfrac>0.)then
            tfrac = 1./sfrac
            if(tsed<0.) ero(gid,:) =  ero(gid,:) + tfrac*sed(:)
            if(tsed>0.) depo(gid,:) =  depo(gid,:) + tfrac*sed(:)
            sumDiff(gid) = sumDiff(gid) + tfrac*tsed
          endif
        elseif(pyZ(gid)<pyZ(ngbid))then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          sfrac = 0.
          sed = 0.
          tsed = 0.
          do r = 1, pyRockNb
            frac = difflay(ngbid,r)/maxlayh(ngbid)
            sfrac = sfrac + frac
            sed(r) = pyCoeff(gid)*frac*flx
            tsed = tsed + sed(r)
          enddo
          if(sfrac>0.)then
            tfrac = 1./sfrac
            if(tsed<0.) ero(gid,:) =  ero(gid,:) + tfrac*sed(:)
            if(tsed>0.) depo(gid,:) =  depo(gid,:) + tfrac*sed(:)
            sumDiff(gid) = sumDiff(gid) + tfrac*tsed
          endif
        endif
      elseif(pyBord(ngbid)<1)then
        if(pyZ(gid)>pyZ(ngbid))then
          flx = pyEdge(gid,p)*(pyZ(ngbid)-pyZ(gid))/pyDist(gid,p)
          sfrac = 0.
          sed = 0.
          tsed = 0.
          do r = 1, pyRockNb
            frac = difflay(gid,r)/maxlayh(gid)
            sfrac = sfrac + frac
            sed(r) = pyCoeff(gid)*frac*flx
            tsed = tsed + sed(r)
          enddo
          if(sfrac>0.)then
            tfrac = 1./sfrac
            if(tsed<0.) ero(gid,:) =  ero(gid,:) + tfrac*sed(:)
            if(tsed>0.) depo(gid,:) =  depo(gid,:) + tfrac*sed(:)
            sumDiff(gid) = sumDiff(gid) + tfrac*tsed
          endif

        endif
      endif
    enddo loop
  endif
  enddo

  return

end subroutine diffsedhillslope

subroutine slumpero(pyStack, pyRcv, pyXY, pyElev, pySfail, borders,pyEro, pylNodesNb, pygNodesNb)

  use classfv
  implicit none

  integer :: pylNodesNb
  integer :: pygNodesNb

  real(kind=8),intent(in) :: pySfail
  integer,dimension(pylNodesNb),intent(in) :: pyStack
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  integer,dimension(pygNodesNb),intent(in) :: borders

  real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev

  real(kind=8),dimension(pygNodesNb),intent(out) :: pyEro

  integer :: n, donor, recvr
  real(kind=8) :: dh, dhmax
  real(kind=8) :: dist, slope

  pyEro = 0.

  do n = pylNodesNb, 1, -1

  donor = pyStack(n) + 1
  recvr = pyRcv(donor) + 1
  dh = (pyElev(donor) - pyElev(recvr))
  dist = sqrt( (pyXY(donor,1)-pyXY(recvr,1))**2.0 + (pyXY(donor,2)-pyXY(recvr,2))**2.0 )
  slope = dh/dist
  dhmax = 0.
  if(slope>pySfail)then
    dhmax = (pySfail-0.001)*dist
  endif
  if(borders(donor) > 0)then
    pyEro(donor) = -dhmax
  endif

  enddo

  return

end subroutine slumpero

subroutine streampower(sedfluxcrit,pyStack, pyRcv, pitID, pitVol1, pitDrain, pyXY, pyArea, pyMaxH, &
pyMaxD, pyDischarge, pyFillH, pyElev, pyRiv, Cero, actlay, perc_dep, slp_cr, sea, db, dt, &
borders, pyDepo, pyEro, sedFluxes, slope, pyDensity, pylNodesNb, pygNodesNb, pyRockNb)

  use classfv
  implicit none

  integer :: pylNodesNb
  integer :: pygNodesNb
  integer :: pyRockNb
  real(kind=8),intent(in) :: dt
  real(kind=8),intent(in) :: sea
  real(kind=8),intent(in) :: db
  real(kind=8),intent(in) :: perc_dep
  real(kind=8),intent(in) :: slp_cr
  real(kind=8),intent(in) :: sedfluxcrit
  integer,dimension(pylNodesNb),intent(in) :: pyStack
  integer,dimension(pygNodesNb),intent(in) :: pyRcv
  integer,dimension(pygNodesNb),intent(in) :: pitID
  integer,dimension(pygNodesNb),intent(in) :: borders
  integer,dimension(pygNodesNb),intent(in) :: pitDrain
  real(kind=8),dimension(pygNodesNb,2),intent(in) :: pyXY
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyArea
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyDischarge
  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(in) :: Cero
  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(in) :: actlay
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyMaxH
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyMaxD
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyFillH
  real(kind=8),dimension(pygNodesNb),intent(in) :: pyElev
  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(in) :: pyRiv
  real(kind=8),dimension(pygNodesNb),intent(in) :: pitVol1

  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(out) :: pyDepo
  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(out) :: pyEro
  real(kind=8),dimension(pygNodesNb,pyRockNb),intent(out) :: sedFluxes
  real(kind=8),dimension(pygNodesNb),intent(out) :: slope
  real(kind=8),dimension(pygNodesNb),intent(out) :: pyDensity

  integer :: n, donor, recvr, nID, tmpID, r
  real(kind=8) :: maxh, dh, waterH, fct, Qt, totflx, totspl, newdist,rhosed,rhowat,tauratio
  real(kind=8) :: dist, slp, slpdh, updh, tmpdist, totdist, width, frac, upperslp, bedfrac
  real(kind=8),dimension(pyRockNb) :: SPL, Qs, Qb, frck, erodep, pitDep
  real(kind=8),dimension(pygNodesNb) :: upZ, updist, pitVol,hypyc
  real(kind=8),dimension(pygNodesNb,pyRockNb) :: bedFluxes

  pyDepo = 0.
  pyEro = 0.
  pyDensity = 1000.
  rhosed = 2000.
  rhowat = 1000.
  hypyc = 0.
  slope = 0.
  pitVol = pitVol1
  sedFluxes = pyRiv * dt
  if(bedslptype > 0) bedFluxes = 0.
  upZ = 1.e6
  updist = 0.

  do n = pylNodesNb, 1, -1

    SPL = 0.
    donor = pyStack(n) + 1
    recvr = pyRcv(donor) + 1
    dh = 0.95*(pyElev(donor) - pyElev(recvr))

    if(pyElev(donor) > sea .and. pyElev(recvr) < sea) dh = 0.99*(pyElev(donor) - sea)
    if( dh < 0.001 ) dh = 0.
    waterH = pyFillH(donor)-pyElev(donor)
    dist = sqrt( (pyXY(donor,1)-pyXY(recvr,1))**2.0 + (pyXY(donor,2)-pyXY(recvr,2))**2.0 )

    ! Compute stream power law
    slpdh = 0.
    bedfrac = 0.02
    totspl = 0
    totdist = 0.
    if( recvr /= donor .and. dh > 0.)then
      ! In case where there is no depression
      if(waterH == 0. )then
        totflx=0.
        do r=1, pyRockNb
          totflx=totflx+sedFluxes(donor,r)
        enddo
        pyDensity(donor) = (totflx/(dt*pyDischarge(donor)))*rhosed+(1-totflx/(dt*pyDischarge(donor)))*rhowat
        if (pyElev(donor)<=sea ) then
           if (pyDensity(donor) >= sedfluxcrit) then
             hypyc(donor) = 1.
             hypyc(recvr) = 1.
           endif
        endif
        if (pyElev(donor) <= sea .and. hypyc(donor) >=.1 ) then
           tauratio = (rhosed-rhowat)/rhosed
        elseif(pyElev(donor) <= sea .and. hypyc(donor)<=0.) then
           tauratio = 0.
        else
           tauratio = 1.
        endif
        if( pyElev(donor)<db) tauratio = 0.

        slope(donor) = dh/dist
        slp = slope(donor)*tauratio
        ! Check if this is an alluvial plain in which case we force deposition
        if(updist(donor) > 0. .and. dist > 0. .and. slp_cr > 0.)then
          updh = upZ(donor) - pyElev(donor)
          if(maxval(sedFluxes(donor,:)) > 0. .and. updh/updist(donor) < slp_cr .and. slp < slp_cr .and. updh > 0)then
            slpdh = perc_dep * updh
            slpdh = min(slpdh,pyMaxD(donor))
          endif
        elseif(incisiontype > 0 .and. dist > 0. .and. updist(donor) > 0.)then
          slpdh = upZ(donor) - pyElev(donor)
        endif

        if(bedslptype > 0 .and. updist(donor) > 0.)then
          ! Compute upper slope
          upperslp = abs(upZ(donor) - pyElev(donor))/updist(donor)
          ! Find bedload fraction in current node
          if(upperslp >= 1./sqrt(3.))then
            bedfrac = 1.
          elseif(bedslptype == 1)then
            bedfrac = 0.98 * sqrt(3.) * upperslp + 0.02
          elseif(bedslptype == 2)then
            bedfrac = (1. / (1. + abs((upperslp - 0.60965)/0.08)**(1.912)) - 0.0201)*1.181 + 0.02
          elseif(bedslptype == 3)then
            bedfrac = (0.8499389 - 1./(1. + abs((upperslp+0.0323)/0.08)**(1.912)))*1.181 + 0.02
          endif
        elseif(bedslptype > 0 .and. updist(donor) == 0.)then
          bedfrac = 0.02
        elseif(bedslptype == 0)then
          bedfrac = 1.
        endif

        ! Compute the stream power law expressed in m/y
        if(dist > 0.)then
          SPL = 0.
          totspl = 0

          ! Get fraction of each rock type present in the active layer
          totflx = 0.
          if(pyRockNb>1)then
            do r = 1, pyRockNb
              totflx = totflx+actlay(donor,r)
            enddo
            frck(1:pyRockNb) = actlay(donor,1:pyRockNb)/totflx
          else
            frck(1) = 1.
          endif

          ! Incision rule types
          ! Detachment limited
          if(incisiontype==0 .and. slpdh == 0.)then
            do r = 1, pyRockNb
              SPL(r) = -Cero(donor,r) * frck(r) * bedfrac * (pyDischarge(donor))**spl_m * (slp)**spl_n
              totspl = totspl + SPL(r)
            enddo
            if(-totspl*dt>dh)then
              if(dh==0.)then
                SPL = 0.
                totspl = 0.
              else
                frac = dh/(-totspl*dt)
                SPL = SPL*frac
                totspl = -dh
              endif
            endif
            if(pyElev(donor)<db)then
              SPL = 0.
              totspl = 0.
            endif

          ! Generalised undercapacity model (linear sedflux dependency)
          elseif(incisiontype==1)then
            Qt = sed_kt * (pyDischarge(donor))**sed_mt * (slp)**sed_nt
            totflx = 0.
            do r = 1, pyRockNb
              if(bedslptype > 0)then
                totflx = totflx + bedFluxes(donor,r)
              else
                totflx = totflx + sedFluxes(donor,r)
              endif
            enddo
            if(Qt>0.)then
              fct = 1. - totflx/Qt
              if(fct<0.) fct = 0.
              if(fct>1.) fct = 1.
            else
              fct = 0.
            endif
            do r = 1, pyRockNb
              SPL(r) = -Cero(donor,r) * frck(r) * fct * (pyDischarge(donor))**spl_m * (slp)**spl_n
              totspl = totspl + SPL(r)
            enddo
            if(-totspl*dt>dh)then
              if(dh==0.)then
                SPL = 0.
                totspl = 0.
              else
                frac = dh/(-totspl*dt)
                SPL = SPL*frac
                totspl = -dh
              endif
            endif

          ! Almost parabolic sedflux dependency
          elseif(incisiontype==2)then
            Qt = sed_kt * (pyDischarge(donor))**sed_mt * (slp)**sed_nt
            totflx = 0.
            do r = 1, pyRockNb
              if(bedslptype > 0)then
                totflx = totflx + bedFluxes(donor,r)
              else
                totflx = totflx + sedFluxes(donor,r)
              endif
            enddo
            if(Qt>0.)then
              frac = totflx/Qt
              if(frac<0.1)then
                fct = 2.6*frac + 0.1
              else
                fct = 1. - 4*(frac-0.5)**2.
              endif
              if(fct<0.) fct = 0.
              if(fct>1.) fct = 1.
            else
              fct = 0.
            endif
            do r = 1, pyRockNb
              SPL(r) = -Cero(donor,r) * frck(r) * fct * (pyDischarge(donor))**spl_m * (slp)**spl_n
              totspl = totspl + SPL(r)
            enddo
            if(-totspl*dt>dh)then
              if(dh==0.)then
                SPL = 0.
                totspl = 0.
              else
                frac = dh/(-totspl*dt)
                SPL = SPL*frac
                totspl = -dh
              endif
            endif

          ! Almost parabolic sedflux dependency
          elseif(incisiontype==3)then
            Qt = sed_kt * (pyDischarge(donor))**sed_mt * (slp)**sed_nt
            totflx = 0.
            do r = 1, pyRockNb
              if(bedslptype > 0)then
                totflx = totflx + bedFluxes(donor,r)
              else
                totflx = totflx + sedFluxes(donor,r)
              endif
            enddo
            if(Qt>0.)then
              frac = totflx/Qt
              if(frac<0.35)then
                fct = exp(-(frac - 0.35)**2/(0.22)**2)
              else
                fct = exp(-(frac - 0.35)**2/(0.6)**2)
              endif
              if(fct<0.) fct = 0.
              if(fct>1.) fct = 1.
            else
              fct = 0.
            endif
            do r = 1, pyRockNb
              SPL(r) = -Cero(donor,r) * frck(r) * fct * (pyDischarge(donor))**spl_m * (slp)**spl_n
              totspl = totspl + SPL(r)
            enddo
            if(-totspl*dt>dh)then
              if(dh==0.)then
                SPL = 0.
                totspl = 0.
              else
                frac = dh/(-totspl*dt)
                SPL = SPL*frac
                totspl = -dh
              endif
            endif

          ! Saltation abrasion incision model
          elseif(incisiontype==4)then
            Qt = sed_kt * (pyDischarge(donor))**sed_mt * (slp)**sed_nt
            totflx = 0.
            do r = 1, pyRockNb
              if(bedslptype > 0)then
                totflx = totflx + bedFluxes(donor,r)
              else
                totflx = totflx + sedFluxes(donor,r)
              endif
            enddo
            if(Qt>0.)then
              fct = 1. - totflx/Qt
              if(fct<0.) fct = 0.
              if(fct>1.) fct = 1.
            else
              fct = 0.
            endif
            ! Channel width
            width = width_kw * (pyDischarge(donor))**width_b
            if(width>0)then
              do r = 1, pyRockNb
                SPL(r) = -Cero(donor,r) * frck(r) * totflx/(dt * width) * fct * (pyDischarge(donor))**spl_m * (slp)**spl_n
                totspl = totspl + SPL(r)
              enddo
            else
              SPL = 0.
            endif
          endif
        endif
      endif
    endif

    maxh = pyMaxH(donor)
    if(waterH > 0.)then
      maxh = waterH
    elseif(pyElev(donor) < sea)then
      maxh = sea - pyElev(donor)
    elseif(slpdh > 0. .and. slp_cr > 0.)then
      maxh = slpdh
    elseif(slpdh > 0. .and. incisiontype > 0)then
      maxh = slpdh
    endif
    maxh = 0.95*maxh

    Qs = 0.
    if(bedslptype > 0) Qb = 0.
    erodep = 0.
    pitDep = 0.
    ! Erosion case
    if(totspl < 0.)then
      ! Sediment volume [m3]
      ! Limit erosion based on active layer rock proportion
      if(pyRockNb>1)then
        do r = 1, pyRockNb
          if(-SPL(r)*dt>actlay(donor,r))then
            erodep(r) = -actlay(donor,r) * pyArea(donor)
          else
            erodep(r) = SPL(r) * dt * pyArea(donor)
          endif
          Qs(r) = -erodep(r) + sedFluxes(donor,r)
          if(bedslptype > 0) Qb(r) = -erodep(r)*bedfrac + bedFluxes(donor,r)
        enddo
      else
        erodep(1) = SPL(1) * dt * pyArea(donor)
        Qs(1) = -erodep(1) + sedFluxes(donor,1)
        if(bedslptype > 0) Qb(1) = -erodep(1)*bedfrac + bedFluxes(donor,1)
      endif

    ! Deposition case
    elseif( totspl >= 0. .and. pyArea(donor) > 0.)then
      ! Fill depression
      if(waterH > 0. )then
        Qs = 0.
        if(bedslptype > 0) Qb = 0.
        erodep = 0.
        totdist = 0.
        do r = 1, pyRockNb
          pitDep(r) = sedFluxes(donor,r)
          totdist = totdist + pitDep(r)
        enddo

      ! Marine deposit
      ! if hypopycnal flow
      elseif(pyElev(donor) <= sea .and. hypyc(donor) <= 0.)then
        ! Add all sediment to the node
        do r = 1, pyRockNb
          erodep(r) = sedFluxes(donor,r)
        enddo
        Qs = 0.
        if(bedslptype > 0) Qb = 0.
      ! Alluvial plain deposit
      !elseif(maxh > 0. .and. waterH == 0. .and. donor /= recvr .and. pyElev(donor) > sea)then
      elseif(maxh > 0. .and. waterH == 0. .and. donor /= recvr )then
        totflx = 0.
        do r = 1, pyRockNb
          totflx = totflx+sedFluxes(donor,r)
        enddo
        if(totflx/pyArea(donor) < maxh)then
          do r = 1, pyRockNb
            erodep(r) = sedFluxes(donor,r)
          enddo
          Qs = 0.
          if(bedslptype > 0) Qb = 0.
        else
          do r = 1, pyRockNb
            frac =  sedFluxes(donor,r)/totflx
            erodep(r) = frac*maxh*pyArea(donor)
            Qs(r) = sedFluxes(donor,r) - erodep(r)
            if(bedslptype > 0) Qb(r) = max(0.,bedFluxes(donor,r) - erodep(r))
          enddo
        endif

      ! Base-level (sink)
      elseif(donor == recvr .and. pyArea(donor) > 0.)then
        do r = 1, pyRockNb
          erodep(r) = sedFluxes(donor,r)
        enddo
        Qs = 0.
        if(bedslptype > 0) Qb = 0.
      else
        erodep = 0.
        do r = 1, pyRockNb
          Qs(r) = sedFluxes(donor,r)
          if(bedslptype > 0) Qb(r) = bedFluxes(donor,r)
        enddo
      endif
    endif

    ! Update sediment volume in receiver node
    if(maxval(pitDep)==0.)then
      do r = 1, pyRockNb
        sedFluxes(recvr,r) = sedFluxes(recvr,r) + Qs(r)
        if(bedslptype > 0) bedFluxes(recvr,r) = bedFluxes(recvr,r) + Qb(r)
        if(erodep(r)<0.)then
          pyEro(donor,r) = pyEro(donor,r) + erodep(r)
        else
          pyDepo(donor,r) = pyDepo(donor,r) + erodep(r)
        endif
      enddo

    ! In case we fill a depression
    elseif(maxval(pitDep)>0. .and. pyArea(pitID(donor)+1)>0.)then
      ! Perform distribution
      tmpID = pitID(donor) + 1

      do while(totdist > 0.)
        ! Get the volume already deposited on the considered node
        tmpdist = 0.
        do r = 1, pyRockNb
          tmpdist = tmpdist + pyDepo(tmpID,r)
        enddo

        ! In case the depression is underwater
        if(pyfillH(tmpID)<sea)then
          if(pyElev(donor)<sea)then
            do r = 1, pyRockNb
              pyDepo(donor,r) = pyDepo(donor,r) + pitDep(r)
            enddo
            totdist = 0.
          else
            do r = 1, pyRockNb
              sedFluxes(recvr,r) = sedFluxes(recvr,r) + pitDep(r)
            enddo
            totdist = 0.
          endif
          nID = recvr

        ! In case the depression is not filled
        elseif(tmpdist+totdist<=pitVol(tmpID))then
          do r = 1, pyRockNb
            pyDepo(tmpID,r) = pyDepo(tmpID,r) + pitDep(r)
          enddo
          totdist = 0.
          nID = tmpID

        ! In case this is an internally drained depression
        elseif(pitDrain(tmpID)+1==tmpID)then
          do r = 1, pyRockNb
            pyDepo(tmpID,r) = pyDepo(tmpID,r) + pitDep(r)
          enddo
          totdist = 0.
          nID = tmpID

        ! Otherwise get the amount to distibute towards draining basins
        else
          if(borders(tmpID) == 0)then
             totdist = 0.
             nID = tmpID
          elseif(tmpdist == pitVol(tmpID))then
             nID = tmpID
          else
             newdist = 0.
             totflx = 0.
             do r = 1, pyRockNb
               frac = pitDep(r)/totdist
               pyDepo(tmpID,r) = pyDepo(tmpID,r) + (pitVol(tmpID) - tmpdist)*frac
               pitDep(r) = (totdist - (pitVol(tmpID) - tmpdist))*frac
               newdist = newdist + pitDep(r)
               totflx = totflx + pyDepo(tmpID,r)
             enddo
             totdist = newdist
             pitVol(tmpID) = totflx
             nID = tmpID
          endif
        endif
        tmpID = pitDrain(nID) + 1
      enddo
    endif

    ! For alluvial deposition
    upZ(recvr) = min(pyElev(donor),upZ(recvr))
    if(upZ(recvr)==pyElev(donor)) updist(recvr) = dist

  enddo

  return

end subroutine streampower


subroutine getid1(volc,vol,alldrain,pit,sumvol,ids,ids2,newNb,newNb2,ptsNb,sedNb)

  use classfv
  implicit none

  integer :: ptsNb
  integer :: sedNb
  integer,dimension(ptsNb),intent(in) :: pit
  integer,dimension(ptsNb),intent(in) :: alldrain
  real(kind=8),dimension(ptsNb,sedNb),intent(in) :: volc
  real(kind=8),dimension(ptsNb),intent(in) :: vol

  integer, intent(out) :: newNb
  integer, intent(out) :: newNb2
  integer,dimension(ptsNb), intent(out) :: ids
  integer,dimension(ptsNb), intent(out) :: ids2
  real(kind=8),dimension(ptsNb), intent(out) :: sumvol

  integer :: p,s

  newNb = 0
  newNb2 = 0
  sumvol = 0.
  ids = 0
  ids2 = 0

  do p = 1, ptsNb
  do s = 1, sedNb
    sumvol(p) = sumvol(p)+volc(p,s)
  enddo
  if(sumvol(p)>vol(p) .and. vol(p)>0.)then
    newNb = newNb+1
    ids(newNb) = p-1
  endif
  if(pit(p)>=0.and.alldrain(p)==pit(p))then
    newNb2 = newNb2+1
    ids2(newNb2) = p-1
  endif
  enddo

  return

end subroutine getid1


subroutine getids(fillH,elev,depo,vol,seal,ids,ids2,ids3,perc,newNb,newNb2,newNb3,ndepo,ptsNb,sedNb)

  use classfv
  implicit none

  integer :: ptsNb
  integer :: sedNb

  real(kind=8),intent(in) :: seal
  real(kind=8),dimension(ptsNb),intent(in) :: fillH
  real(kind=8),dimension(ptsNb),intent(in) :: elev
  real(kind=8),dimension(ptsNb),intent(in) :: vol
  real(kind=8),dimension(ptsNb,sedNb),intent(in) :: depo

  integer, intent(out) :: newNb
  integer, intent(out) :: newNb2
  integer, intent(out) :: newNb3
  integer,dimension(ptsNb), intent(out) :: ids
  integer,dimension(ptsNb), intent(out) :: ids2
  integer,dimension(ptsNb), intent(out) :: ids3
  real(kind=8),dimension(ptsNb,sedNb), intent(out) :: perc
  real(kind=8),dimension(ptsNb,sedNb), intent(out) :: ndepo

  integer :: p,s,in
  real(kind=8) :: sumdep, sumperc, nvol

  newNb = 0
  newNb2 = 0
  newNb3 = 0
  ids = 0
  ids2 = 0
  ids3 = 0
  perc = 0.
  ndepo = depo

  do p = 1, ptsNb
  in = 0
  sumdep = 0.
  sumperc = 0.
  ! Get alluvial plain deposition ID
  if(elev(p)>seal .and. fillH(p)==elev(p))then
    in = 1
    do s = 1, sedNb
      sumdep = sumdep+depo(p,s)
    enddo
    if(sumdep>0.)then
      newNb = newNb+1
      ids(newNb) = p-1
    endif
  endif

  ! Get land pit deposition ID
  if(elev(p)>seal .and. fillH(p)>seal .and. vol(p)>0.)then
    if(sumdep==0. .and. in == 0)then

      nvol = 0.
      do s = 1, sedNb
        if(depo(p,s)/vol(p)>0.)then
          nvol = nvol+depo(p,s)
        endif
      enddo
      do s = 1, sedNb
        if(depo(p,s)/vol(p)>0.)then
          perc(p,s) = depo(p,s)/nvol
          sumperc = sumperc+perc(p,s)
        else
          ndepo(p,s) = 0.
          perc(p,s) = 0.
        endif
      enddo
    endif
    if(sumperc>0.)then
      newNb2 = newNb2+1
      ids2(newNb2) = p-1
      if(abs(sumperc-1.)>1.e-8)then
        perc(p,:) = perc(p,:)/sumperc
      endif
    endif
  endif

  ! Get water deposition ID
  if(elev(p)<=seal)then
    do s = 1, sedNb
      sumdep = sumdep+depo(p,s)
    enddo
    if(sumdep>0.)then
      newNb3 = newNb3+1
      ids3(newNb3) = p-1
    endif
  endif
  enddo

  return

end subroutine getids
