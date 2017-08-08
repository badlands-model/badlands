!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm
module pdstack

  use pdclass

  implicit none

contains

  subroutine initialisePD(elevation, demH, pydnodes)

    integer :: pydnodes, k
    real(kind=8),intent(in) :: elevation(pydnodes)
    real(kind=8),intent(inout) :: demH(pydnodes)

    s1 = 0
    s2 = 0

    demH(1:bds) = elevation(1:bds)
    do k = bds+1,pydnodes
        demH(k) = 1.e6
        s1 = s1 + 1
        data1(s1) = k
    enddo

  end subroutine initialisePD

  subroutine fillPD(elevation, sealevel, demH, pydnodes)

    logical :: change
    integer :: pydnodes, n, k, p
    real(kind=8) :: hmin
    real(kind=8),intent(in) :: sealevel
    real(kind=8),intent(in) :: elevation(pydnodes)
    real(kind=8),intent(inout) :: demH(pydnodes)

    change = .true.
    do while(change)
      change = .false.
      do n = 1, s1
        k = data1(n)
        if( demH(k) > elevation(k) )then
          ! Get minimum value
          hmin = 2.e6
          loop: do p = 1, 20
            if( neighbours(k,p) < 0 ) exit loop
            hmin = min(hmin,demH(neighbours(k,p)+1))
          enddo loop
          if( elevation(k) >= hmin + eps )then
            demH(k) = elevation(k)
          else
            if( demH(k) > hmin + eps )then
              demH(k) = hmin + eps
              if(elevation(k)>=sealevel)then
                if( demH(k) - elevation(k) > fill_TH)then
                  demH(k) = elevation(k) + fill_TH
                else
                  change = .true.
                endif
              else
                if( demH(k) - sealevel > fill_TH)then
                  demH(k) = sealevel + fill_TH
                else
                  change = .true.
                endif
              endif
            endif
            s2 = s2 + 1
            data2(s2) = k
          endif
        endif
      enddo
      if( s2 > 0 )then
        s1 = s2
        data1(1:s1) = data2(1:s2)
        s2 = 0
      endif
    enddo

    return

  end subroutine fillPD

  subroutine allfillPD(elevation, demH, pydnodes)

    logical :: change
    integer :: pydnodes, n, k, p
    real(kind=8) :: hmin
    real(kind=8),intent(in) :: elevation(pydnodes)
    real(kind=8),intent(inout) :: demH(pydnodes)

    change = .true.
    do while(change)
      change = .false.
      do n = 1, s1
        k = data1(n)
        if( demH(k) > elevation(k) )then
          ! Get minimum value
          hmin = 2.e6
          loop: do p = 1, 20
            if( neighbours(k,p) < 0 ) exit loop
            hmin = min(hmin,demH(neighbours(k,p)+1))
          enddo loop
          if( elevation(k) >= hmin + eps )then
            demH(k) = elevation(k)
          else
            if( demH(k) > hmin + eps )then
              demH(k) = hmin + eps
              change = .true.
            endif
            s2 = s2 + 1
            data2(s2) = k
          endif
        endif
      enddo
      if( s2 > 0 )then
        s1 = s2
        data1(1:s1) = data2(1:s2)
        s2 = 0
      endif
    enddo

    return

  end subroutine allfillPD

  subroutine marine_distribution(elevation, seavol, sealevel, border, depIDs, diffsed, pydnodes, pyIDs, pyRockNb)

    integer :: pydnodes, pyIDs, pyRockNb
    integer,dimension(pyIDs),intent(in) :: depIDs
    integer,dimension(pydnodes),intent(in) :: border
    real(kind=8),intent(in) :: sealevel
    real(kind=8),dimension(pydnodes,pyRockNb),intent(in) :: seavol
    real(kind=8),dimension(pydnodes),intent(in) :: elevation

    real(kind=8),dimension(pydnodes,pyRockNb),intent(out) :: diffsed

    real(kind=8),dimension(pydnodes) :: elev, seadep, newelev

    integer :: it, s, m, n, p, k, pid, nid, id, nup, ndown, ngbh(20)
    real(kind=8) :: dh, vol, minz, maxz

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
              vol = max(0.,diffprop*(maxz-elev(id))*area(id))

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

  subroutine pitparams(pyNgbs,pyArea,pyDiff,pyProp,fillTH,epsilon,pybounds,pydnodes)

    integer :: pydnodes
    integer,intent(in) :: pybounds
    real(kind=8),intent(in) :: fillTH
    real(kind=8),intent(in) :: epsilon
    integer,intent(in) :: pyDiff
    integer,intent(in) :: pyNgbs(pydnodes,20)
    real(kind=8),intent(in) :: pyProp
    real(kind=8),intent(in) :: pyArea(pydnodes)

    dnodes = pydnodes

    diffnbmax = pyDiff
    diffprop = pyProp
    bds = pybounds
    block_size = pydnodes - bds
    eps = epsilon
    fill_TH = fillTH

    call defineparameters

    neighbours = pyNgbs
    area = pyArea

    return

  end subroutine pitparams

  subroutine pitfilling(elevation,allfill,sealevel,demH,pydnodes)

    integer :: pydnodes
    integer,intent(in) :: allfill
    real(kind=8),intent(in) :: sealevel
    real(kind=8),intent(in) :: elevation(pydnodes)

    real(kind=8),intent(out) :: demH(pydnodes)

    ! Initialisation phase
    call initialisePD(elevation,demH,pydnodes)

    ! Filling phase
    if(allfill == 0)then
      call fillPD(elevation,sealevel,demH,pydnodes)
    else
      call allfillPD(elevation,demH,pydnodes)
    endif

    return

  end subroutine pitfilling

  subroutine getactlay(alay,layTH,laySD,alayS,nbPts,nbLay,nbSed)

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

  subroutine updatestrati(layS,layH,eros,depo,newH,newS,nbPts,nbLay,nbSed)

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
              if(newS(n,k,s)>=ero)then
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

  subroutine stratcarb(layS,layH,clastic,carb,carb2,pel,newH,newS,nbPts,nbLay,nbSed)

    integer :: nbPts
    integer :: nbLay
    integer :: nbSed
    real(kind=8),dimension(nbPts),intent(in) :: pel
    real(kind=8),dimension(nbPts),intent(in) :: carb
    real(kind=8),dimension(nbPts),intent(in) :: carb2
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
            do s = 1, 3
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
            newS(n,k,1:3) = 0.
          endif
        enddo lp
      endif
      if(clastic(n)>0.)then
        newH(n,nbLay) = newH(n,nbLay) + clastic(n)
        newS(n,nbLay,1) = newS(n,nbLay,1) + clastic(n)
      endif
      if(carb(n)>0.)then
        newH(n,nbLay) = newH(n,nbLay) + carb(n)
        newS(n,nbLay,2) = newS(n,nbLay,2) + carb(n)
      endif
      if(carb2(n)>0.)then
        newH(n,nbLay) = newH(n,nbLay) + carb2(n)
        newS(n,nbLay,3) = newS(n,nbLay,3) + carb2(n)
      endif
      if(pel(n)>0.)then
        newH(n,nbLay) = newH(n,nbLay) + pel(n)
        newS(n,nbLay,4) = newS(n,nbLay,4) + pel(n)
      endif
    enddo

    return

  end subroutine stratcarb

end module pdstack
