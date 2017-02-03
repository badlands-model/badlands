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

  subroutine marine_sed(elevation, seavol, border, sealevel, seadep, pydnodes)

    integer :: pydnodes
    integer,dimension(pydnodes),intent(in) :: border

    real(kind=8),intent(in) :: sealevel
    real(kind=8),dimension(pydnodes),intent(in) :: elevation
    real(kind=8),dimension(pydnodes),intent(in) :: seavol

    real(kind=8),dimension(pydnodes),intent(out) :: seadep


    real(kind=8),dimension(pydnodes) :: difo, difp, cdif, elev, depo

    integer :: n, k, p, iter, ndown
    integer,dimension(20) :: down

    real(kind=8) :: cc, volsum, difsed, fc, difmax
    real(kind=8),dimension(20) :: vol

    seadep = 0.
    elev = elevation

    depo = 0.
    ! Deposit sediment in diff_nb increments
    do n = 1, diff_nb

      ! Initialise sediments
      difsed = 0.
      difo = seavol/diff_nb
      difp = 0.

      ! Move sediment to downstream nodes as long as there is sediment in the buffer.
      difmax = 1.0e20
      iter = 0
      sediment_diffusion: do
        if(difmax<=diff_res*meanarea .or. iter>=max_it_cyc) exit sediment_diffusion
        difsed = 0.
        difmax = 0.
        iter = iter + 1
        ! Determine the fraction of buffer that will be deposited at current node
        cdif = 0.
        ! Get maximum elevation at current nodes to ensure sediment stability
        call find_maximum_elevation(difo,elev,sealevel,border,cdif,depo)

        ! Calculate individual and cumulative accomodation space for each downstream node
        do k = 1,pydnodes
          if(difo(k)>0.)then
            if(cdif(k)<1.)then
              call get_available_volume_remaining(k,elev,volsum,vol,down,ndown)
              ! Remaining volume of sediment to distribute
              cc = difo(k)*(1.-cdif(k))
              difsed = difsed + cc

              ! If enough downstream volume is available weight sediment by
              ! available volume and diffuse sediment to downstream nodes
              if(volsum>=cc)then
                ! Loop over neighboring cell
                loop1: do p = 1, 20
                  if(neighbours(k,p)<0) exit loop1
                  difp(neighbours(k,p)+1) = difp(neighbours(k,p)+1)+cc*vol(p)/volsum
                enddo loop1

              ! If not enough downstream volume is available weight sediment
              ! by available volume and distribute rest equally among
              ! connected downstream nodes or if no downstream nodes are
              ! connected assume that enough downstream volume is available
              ! and weight sediment by the number of neighboring nodes
              else
                if(ndown>0)then
                  fc = (cc-volsum)/float(ndown)
                  ! Loop over neighboring cell
                  loop2: do p = 1, 20
                    if(neighbours(k,p)<0) exit loop2
                    if(down(p)==1)then
                      difp(neighbours(k,p)+1) = difp(neighbours(k,p)+1)+vol(p)+fc
                    endif
                  enddo loop2
                ! If there are no downtream nodes, just smear it all around
                ! and hope the next iteration looks after it
                else
                  fc = cc/maxngb(k)
                  ! Loop over neighboring cell
                  loop3: do p = 1, 20
                     if(neighbours(k,p)<0) exit loop3
                     difp(neighbours(k,p)+1) = difp(neighbours(k,p)+1)+fc
                  enddo loop3
                endif
              endif
            endif
          endif
        enddo

        ! Store sediment still to be diffused in difo for next iteration
        difo = difp
        difp = 0.0_8
        difmax = max(difsed,difmax)

      enddo sediment_diffusion
    enddo

    seadep = depo

    return

  end subroutine marine_sed

  subroutine pitparams(pyNgbs,pyDist,pyArea,pySlp,fillTH,epsilon,pybounds,pydnodes)

    integer :: pydnodes
    integer,intent(in) :: pybounds
    real(kind=8),intent(in) :: fillTH
    real(kind=8),intent(in) :: epsilon
    real(kind=8),intent(in) :: pySlp
    integer,intent(in) :: pyNgbs(pydnodes,20)
    real(kind=8),intent(in) :: pyDist(pydnodes,20)
    real(kind=8),intent(in) :: pyArea(pydnodes)

    integer :: k, p

    dnodes = pydnodes

    marineslope = pySlp
    bds = pybounds
    block_size = pydnodes - bds
    eps = epsilon
    fill_TH = fillTH

    call defineparameters

    neighbours = pyNgbs
    edge_dist = pyDist
    area = pyArea
    meanarea = 0.

    do k = 1,dnodes
      maxngb(k) = 0
      meanarea = meanarea+area(k)
      loop: do p = 1, 20
        if(neighbours(k,p)<0) exit loop
        maxngb(k) = maxngb(k) + 1
      enddo loop
    enddo
    meanarea = meanarea/dnodes

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
