!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm
subroutine marine_distribution(elevation, seavol, sealevel, border, depIDs, pySlp, diffsed, pydnodes, pyIDs, pyRockNb)

  use classpd
  implicit none

  integer :: pydnodes, pyIDs, pyRockNb
  integer,dimension(pyIDs),intent(in) :: depIDs
  integer,dimension(pydnodes),intent(in) :: border
  real(kind=8),intent(in) :: sealevel
  real(kind=8),dimension(pydnodes,pyRockNb),intent(in) :: seavol
  real(kind=8),dimension(pydnodes),intent(in) :: elevation
  real(kind=8),dimension(pydnodes),intent(in) :: pySlp

  real(kind=8),dimension(pydnodes,pyRockNb),intent(out) :: diffsed

  real(kind=8),dimension(pydnodes) :: elev, seadep, newelev

  integer :: it, s, m, n, p, k, pid, nid, id, nup, ndown, ngbh(20)
  real(kind=8) :: dh, vol, minz, maxz, dprop

  dnodes = pydnodes
  elev = elevation

  do s = 1, pyRockNb
    newelev = elev
    do m = 1, diffnbmax
      seadep = seavol(:,s)/float(diffnbmax)
      do k = 1, pyIDs
        n = depIDs(k)+1
        id = n
        pid = n
        it = 0
        sfd_loop: do
            if(border(id)<1)then
              seadep(id) = 0.
              exit sfd_loop
            endif
            maxz = -1.e8
            loop0: do p = 1, 20
              if( neighbours(id,p) < 0 ) exit loop0
              if(maxz<elev(neighbours(id,p)+1)) maxz = elev(neighbours(id,p)+1)
            enddo loop0
            if(maxz>sealevel) maxz = sealevel
            if(maxz<elev(id)) maxz = elev(id)
            if(propA+propB == 0)then
              dprop = diffprop
            else
              dprop = ((0.9_8)/(1.0_8+exp(propA*(pySlp(id))))+propB)
            endif
            vol = max(0.,dprop*(maxz-elev(id))*area(id))

            if(it>max_it_cyc)then
              elev(id) = elev(id) + seadep(id)/area(id)
              seadep(id) = 0.
              exit sfd_loop
            endif
            if(seadep(id)/area(id)<diff_res)then
              elev(id) = elev(id) + seadep(id)/area(id)
              seadep(id) = 0.
              exit sfd_loop
            endif
            it = it+1
            if(seadep(id)<vol)then
              elev(id) = elev(id) + seadep(id)/area(id)
              seadep(id) = 0.
              exit sfd_loop
            else
              seadep(id) = seadep(id) - vol
              elev(id) = elev(id) + vol/area(id)
            endif
            if(seadep(id)>0.)then
              minz = elev(id)
              nid = 0
              nup = 0
              ndown = 0
              maxz = -1.e8
              loop: do p = 1, 20
                if( neighbours(id,p) < 0 ) exit loop
                ndown = ndown + 1
                ngbh(ndown) = neighbours(id,p)+1
                if(minz>elev(neighbours(id,p)+1))then
                  nid = neighbours(id,p)+1
                  minz = elev(nid)
                endif
                if(maxz<elev(neighbours(id,p)+1))then
                  maxz = elev(neighbours(id,p)+1)
                  nup = neighbours(id,p)+1
                endif
              enddo loop
              if(nid==0)then
                dh = maxz-elev(id)+0.1
                if(seadep(id) > dh*area(id))then
                  elev(id) = elev(id) + dh
                  seadep(id) = seadep(id) - dh*area(id)
                  pid = depIDs(k)+1
                  nid = depIDs(k)+1
                  seadep(nid) = seadep(nid)+seadep(id)
                  seadep(id) = 0.
                  id = nid
                  it = 0
                else
                  elev(id) = elev(id) + seadep(id)/area(id)
                  seadep(id) = 0.
                  exit sfd_loop
                endif
              else
                if(minz>elev(id))then
                  dh = (minz-elev(id))+0.1
                  if(seadep(id) > dh*area(id))then
                    elev(id) = elev(id) + dh
                    seadep(id) = seadep(id) - dh*area(id)
                  else
                    elev(id) = elev(id) + seadep(id)/area(id)
                    seadep(id) = 0.
                    exit sfd_loop
                  endif
                  pid = depIDs(k)+1
                  nid = depIDs(k)+1
                  seadep(nid) = seadep(nid)+seadep(id)
                  seadep(id) = 0.
                  id = nid
                  it = 0
                else
                  seadep(nid) = seadep(nid)+seadep(id)
                  seadep(id) = 0.
                  pid = id
                  id = nid
                endif

              endif
            else
              seadep(id) = 0.
              exit sfd_loop
            endif
        enddo sfd_loop

      enddo
    enddo
    diffsed(:,s) = elev - newelev
  enddo

  return

end subroutine marine_distribution

subroutine pitparams(pyNgbs,pyArea,pyDiff,pyProp,pyPropa,pyPropb,fillTH,epsilon,pybounds,pydnodes)

  use classpd
  implicit none

  integer :: pydnodes
  integer,intent(in) :: pybounds
  real(kind=8),intent(in) :: fillTH
  real(kind=8),intent(in) :: epsilon
  integer,intent(in) :: pyDiff
  integer,intent(in) :: pyNgbs(pydnodes,20)
  real(kind=8),intent(in) :: pyProp
  real(kind=8),intent(in) :: pyPropa
  real(kind=8),intent(in) :: pyPropb
  real(kind=8),intent(in) :: pyArea(pydnodes)

  dnodes = pydnodes

  diffnbmax = pyDiff
  diffprop = pyProp
  propA = pyPropa
  propB = pyPropb
  bds = pybounds
  block_size = pydnodes - bds
  eps = epsilon
  fill_TH = fillTH

  call defineparameters

  neighbours = pyNgbs
  area = pyArea

  print *, "propA: ", propA
  print *, "propB: ", propB

  return

end subroutine pitparams

subroutine pitfilling(elevation,allfill,sealevel,demH,pydnodes)

  use classpd
  implicit none

  integer :: pydnodes
  integer,intent(in) :: allfill
  real(kind=8),intent(in) :: sealevel
  real(kind=8),intent(in) :: elevation(pydnodes)

  real(kind=8),intent(out) :: demH(pydnodes)
  call initialisePD(elevation,demH,pydnodes)

  ! Filling phase
  if(allfill == 0)then
    call fillPD(elevation,sealevel,demH,pydnodes)
    ! call fillBarnes(elevation,sealevel,demH,pydnodes)
  else
    ! call fillBarnes(elevation,sealevel,demH,pydnodes)
    ! Initialisation phase
    call allfillPD(elevation,demH,pydnodes)
  endif

  return

end subroutine pitfilling

subroutine getactlay(alay,layTH,laySD,alayS,nbPts,nbLay,nbSed)

  use classpd
  implicit none

  integer :: nbPts
  integer :: nbLay
  integer :: nbSed
  real(kind=8),dimension(nbPts),intent(in) :: alay
  real(kind=8),dimension(nbPts,nbLay),intent(in) :: layTH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(in) :: laySD
  real(kind=8),dimension(nbPts,nbSed),intent(out) :: alayS

  integer :: n,k,s,lid
  real(kind=8) :: cumh,alayh,prop

  alayS = 0.

  do n = 1, nbPts
    ! Compute cumulative stratal thicknesses
    cumh = layTH(n,nbLay)
    alayh = 0.
    lid = nbLay
    if(cumh<=alay(n))then
      alayh = cumh
      alayS(n,:) = laySD(n,nbLay,:)
      lp: do k = nbLay-1,1,-1
        cumh = cumh+layTH(n,k)
        ! Store stratal thicknesses lower than active layer thickness
        if(alay(n)>=cumh)then
          alayh = cumh
          do s = 1, nbSed
            alayS(n,s) = alayS(n,s)+laySD(n,k,s)
          enddo
        else
          lid = k
          exit lp
        endif
      enddo lp
    endif

    ! Find the proportion of sediment that still needs to be passed to the active layer
    if(layTH(n,lid)>0.)then
      prop = (alay(n)-alayh)/layTH(n,lid)
      prop = min(prop,1.)
      prop = max(prop,0.)
      if(prop>0.)then
        do s = 1, nbSed
          alayS(n,s) = alayS(n,s)+prop*laySD(n,lid,s)
        enddo
      endif
    endif
  enddo

  return

end subroutine getactlay

subroutine getactlay2(alay,layTH,laySD,alayS,nbPts,nbLay,nbSed)

  use classpd
  implicit none

  integer :: nbPts
  integer :: nbLay
  integer :: nbSed
  real(kind=8),intent(in) :: alay
  real(kind=8),dimension(nbPts,nbLay),intent(in) :: layTH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(in) :: laySD
  real(kind=8),dimension(nbPts,nbSed),intent(out) :: alayS

  integer :: n,k,s,lid
  real(kind=8) :: cumh,alayh,prop

  alayS = 0.
  do n = 1, nbPts
    ! Compute cumulative stratal thicknesses
    cumh = layTH(n,nbLay)
    alayh = 0.
    lid = nbLay
    if(cumh<=alay)then
      alayh = cumh
      alayS(n,:) = laySD(n,nbLay,:)
      lp: do k = nbLay-1,1,-1
        cumh = cumh+layTH(n,k)
        ! Store stratal thicknesses lower than active layer thickness
        if(alay>=cumh)then
          alayh = cumh
          do s = 1, nbSed
            alayS(n,s) = alayS(n,s)+laySD(n,k,s)
          enddo
        else
          lid = k
          exit lp
        endif
      enddo lp
    endif
    ! Find the proportion of sediment that still needs to be passed to the active layer
    if(layTH(n,lid)>0.)then
      prop = (alay-alayh)/layTH(n,lid)
      prop = min(prop,1.)
      prop = max(prop,0.)
      if(prop>0.)then
        do s = 1, nbSed
          alayS(n,s) = alayS(n,s)+prop*laySD(n,lid,s)
        enddo
      endif
    endif
  enddo

  return

end subroutine getactlay2

subroutine updatestrati(layS,layH,eros,depo,newH,newS,nbPts,nbLay,nbSed)

  use classpd
  implicit none

  integer :: nbPts
  integer :: nbLay
  integer :: nbSed
  real(kind=8),dimension(nbPts,nbSed),intent(in) :: eros
  real(kind=8),dimension(nbPts,nbSed),intent(in) :: depo
  real(kind=8),dimension(nbPts,nbLay),intent(in) :: layH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(in) :: layS

  real(kind=8),dimension(nbPts,nbLay),intent(out) :: newH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(out) :: newS

  integer :: n,k,s
  real(kind=8) :: ero,dep

  newS = layS
  newH = layH
  do n = 1, nbPts
    do s = 1,nbSed
      ero = -eros(n,s)
      dep = depo(n,s)
      if(ero>0.)then
        lp: do k = nbLay,1,-1
          if(newS(n,k,s)>0.)then
            if(newS(n,k,s)>ero)then
              newS(n,k,s) = newS(n,k,s) - ero
              newH(n,k) = newH(n,k) - ero
              exit lp
            else
              ero = ero - newS(n,k,s)
              newH(n,k) = newH(n,k) - newS(n,k,s)
              newS(n,k,s) = 0.
            endif
          endif
        enddo lp
      endif
      if(dep>0.)then
        newH(n,nbLay) = newH(n,nbLay) + dep
        newS(n,nbLay,s) = newS(n,nbLay,s) + dep
      endif
    enddo
  enddo

  return

end subroutine updatestrati

subroutine updatecstrati(layS,layH,eros,depo,newH,newS,nbPts,nbLay,nbSed)

  use classpd
  implicit none

  integer :: nbPts
  integer :: nbLay
  integer :: nbSed
  real(kind=8),dimension(nbPts),intent(in) :: eros
  real(kind=8),dimension(nbPts),intent(in) :: depo
  real(kind=8),dimension(nbPts,nbLay),intent(in) :: layH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(in) :: layS

  real(kind=8),dimension(nbPts,nbLay),intent(out) :: newH
  real(kind=8),dimension(nbPts,nbLay),intent(out) :: newS

  integer :: n,k
  real(kind=8) :: ero,dep

  newS = layS(:,:,1)
  newH = layH

  do n = 1, nbPts
    ero = -eros(n)
    dep = depo(n)
    if(ero>0.)then
      lp: do k = nbLay,1,-1
        if(newS(n,k)>0.)then
          if(newS(n,k)>=ero)then
            newS(n,k) = newS(n,k) - ero
            newH(n,k) = newH(n,k) - ero
            exit lp
          else
            ero = ero - newS(n,k)
            newH(n,k) = newH(n,k) - newS(n,k)
            newS(n,k) = 0.
          endif
        endif
      enddo lp
    endif
    if(dep>0.)then
      newH(n,nbLay) = newH(n,nbLay) + dep
      newS(n,nbLay) = newS(n,nbLay) + dep
    endif
  enddo

  return

end subroutine updatecstrati

subroutine stratcarb(layS,layH,clastic,newH,newS,nbPts,nbLay,nbSed)

  use classpd
  implicit none

  integer :: nbPts
  integer :: nbLay
  integer :: nbSed
  real(kind=8),dimension(nbPts),intent(in) :: clastic
  real(kind=8),dimension(nbPts,nbLay),intent(in) :: layH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(in) :: layS

  real(kind=8),dimension(nbPts,nbLay),intent(out) :: newH
  real(kind=8),dimension(nbPts,nbLay,nbSed),intent(out) :: newS

  integer :: n,k,s
  real(kind=8) :: ero

  newS = layS
  newH = layH

  do n = 1, nbPts
    if(clastic(n)<0.)then
      ero = -clastic(n)
      lp: do k = nbLay,1,-1
        if(newH(n,k)>ero)then
          do s = 1, nbSed
            if(newS(n,k,s)>=ero)then
              newS(n,k,s) = newS(n,k,s) - ero
              newH(n,k) = newH(n,k) - ero
              exit lp
            else
              ero = ero - newS(n,k,s)
              newH(n,k) = newH(n,k) - newS(n,k,s)
              newS(n,k,s) = 0.
            endif
          enddo
        else
          ero = ero - newH(n,k)
          newH(n,k) = 0.
          newS(n,k,1:nbSed) = 0.
        endif
      enddo lp
    endif
    if(clastic(n)>0.)then
      newH(n,nbLay) = newH(n,nbLay) + clastic(n)
      newS(n,nbLay,1) = newS(n,nbLay,1) + clastic(n)
    endif
  enddo

  return

end subroutine stratcarb
