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

  subroutine marine_distribution(elevation, seavol, sealevel, border, depIDs, diffsed, pydnodes, pyIDs)

    integer :: pydnodes, pyIDs
    integer,dimension(pyIDs),intent(in) :: depIDs
    integer,dimension(pydnodes),intent(in) :: border
    real(kind=8),intent(in) :: sealevel
    real(kind=8),dimension(pydnodes),intent(in) :: seavol
    real(kind=8),dimension(pydnodes),intent(in) :: elevation

    real(kind=8),dimension(pydnodes),intent(out) :: diffsed

    real(kind=8),dimension(pydnodes) :: elev, seadep

    integer :: it, m, n, p, k, pid, nid, id, nup, ndown, ngbh(20)
    real(kind=8) :: dh, vol, minz, maxz

    dnodes = pydnodes
    elev = elevation

    do m = 1, diffnbmax
      seadep = seavol/float(diffnbmax)
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
            vol = max(0.,diffprop*(sealevel-elev(id))*area(id))
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
              minz = 1.e8
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
                  ! vol = (seadep(id) - dh*area(id))/ndown
                  ! do p = 1, ndown
                  !   elev(ngbh(p)) = elev(ngbh(p)) + vol/area(id)
                  ! enddo
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
                  seadep(nid) = seadep(id)
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

    diffsed = elev - elevation

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

end module pdstack
