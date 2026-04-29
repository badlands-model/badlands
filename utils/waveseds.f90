!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Reduced complexity model for wave/sedimentation computation

subroutine airymodel(pdx,pdd,ph0,pdepth,psrc,pinland,pshadow,pc,pl,ptravel,pwaveH,pnumrow,pnumcol)

  implicit none

  integer :: pnumrow,pnumcol
  integer,intent(in) :: pshadow
  integer,dimension(pnumrow,pnumcol),intent(in) :: pinland

  real(kind=8),intent(in) :: pdx
  real(kind=8),intent(in) :: pdd
  real(kind=8),intent(in) :: ph0
  real(kind=8),dimension(pnumrow,pnumcol),intent(in) :: pdepth
  real(kind=8),dimension(pnumrow,pnumcol),intent(in) :: psrc

  real(kind=8),dimension(pnumrow,pnumcol),intent(out) :: pc
  real(kind=8),dimension(pnumrow,pnumcol),intent(out) ::  pl
  real(kind=8),dimension(pnumrow,pnumcol),intent(out) ::  ptravel
  real(kind=8),dimension(pnumrow,pnumcol),intent(out) ::  pwaveH

  integer :: i, j, keeploop, ix, jx, k

  real(kind=8) :: TM, MN, l0, f0, k0, dt, tperiod0
  real(kind=8) :: cg0, c0, kh, tmp, frac, n
  real(kind=8),dimension(pnumrow,pnumcol) :: ks

  real(kind=8), parameter::pi    = 3.1415926535897931_8
  real(kind=8), parameter::onpi  = 1.0_8/3.1415926535897931_8
  real(kind=8), parameter::pi2   = pi*2.0_8
  real(kind=8), parameter::pion2 = pi*0.5_8
  real(kind=8), parameter::onpi2 = 1.0_8/pi2
  real(kind=8), parameter::grav  = 9.81_8

  integer,dimension(20)::iradius=(/-2,-1,1,2,-2,-1,0,1,2,-1,0,1,-2,-1,0,1,2,-1,0,-1/)
  integer,dimension(20)::jradius=(/ 0,0,0,0,1,1,1,1,1, 2,2,2,-1,-1,-1,-1,-1,-2,-2,-2/)

  real(kind=8),dimension(20)::dist=(/2.,1.,1.,2.,sqrt(5.),sqrt(2.),1., &
      sqrt(2.),sqrt(5.),sqrt(5.),2.,sqrt(5.),sqrt(5.),sqrt(2.), &
      1.,sqrt(2.),sqrt(5.),sqrt(5.),2.,sqrt(5.)/)

  pc = 0.
  ks = 1.
  pwaveH = 0.

  ! calculate wave length (deep water)
  tperiod0 = max(0.47*ph0+6.76, pi*pi*sqrt(ph0/grav))

  l0=grav*tperiod0**2*onpi2
  f0=pi2/tperiod0
  k0=pi2/l0
  ! airy wave theory, deep water phase speed
  cg0=0.5_8*sqrt(grav/k0)  ! group speed
  c0=grav*tperiod0*onpi2

  ! set the step size
  pl=l0
  do j = 1, pnumcol
    do i = 1, pnumrow
      ! conct contains all areas not in the shadow of land
      ! if areas are "exposed and in deep water, give them the "open water" conditions
      TM=l0
      if(pinland(i,j)==0)then
        do
          MN=0.5_8*(pl(i,j)+TM)
          TM=pl(i,j)
          pl(i,j)=l0*tanh(pi2*pdepth(i,j)/MN)
          if(abs(pl(i,j)-TM)<1.0e-8_8)exit
        enddo
        pc(i,j)=c0*pl(i,j)/l0
        kh = pdepth(i,j)*pi2/pl(i,j)
        tmp = 1.+2.*kh/sinh(2.*kh)
        pwaveH(i,j) = ph0/sqrt(tanh(kh)*tmp)
        n = 0.5*tmp
        ks(i,j) = sqrt(c0/(2.*n*pc(i,j)))
      endif
    end do
  end do

  ! Assign source points
  ptravel = psrc

  ! Perform Huygen's principle to find travel time and wave front
  keeploop = 1
  do
      keeploop = 0
      do j = 1, pnumcol
        do i = 1, pnumrow
            if(ptravel(i,j)>=0)then
                do k = 1, 20
                    ix = i+iradius(k)
                    jx = j+jradius(k)
                    if(ix>0 .and. ix<=pnumrow .and. jx>0 .and. jx<=pnumcol)then
                        if(pinland(ix,jx)==1)then
                            ptravel(ix,jx)=-1
                        else if(ptravel(ix,jx)<0)then
                            ptravel(ix,jx) = ptravel(i,j)+dist(k)*pdx/pc(i,j)
                            keeploop = 1
                            if(pdepth(i,j)/pl(i,j)<0.5)then
                                frac = 2.*(1.-pdd)*pdepth(i,j)/pl(i,j)+pdd
                                if(pwaveH(ix,jx)>frac*pwaveH(i,j)) pwaveH(ix,jx)=frac*pwaveH(i,j)
                            else
                                if(pshadow==1 .and. pwaveH(ix,jx)>pwaveH(i,j)) pwaveH(ix,jx)=pwaveH(i,j)
                            endif
                        else
                            dt = ptravel(i,j)+dist(k)*pdx/pc(i,j)
                            if(ptravel(ix,jx)>dt .and. dt>0)then
                                ptravel(ix,jx) = dt
                                if(pdepth(i,j)/pl(i,j)<0.5)then
                                    frac = 2.*(1.-pdd)*pdepth(i,j)/pl(i,j)+pdd
                                    if(pwaveH(ix,jx)>frac*pwaveH(i,j)) pwaveH(ix,jx)=frac*pwaveH(i,j)
                                else
                                    if(pshadow==1 .and. pwaveH(ix,jx)>pwaveH(i,j)) pwaveH(ix,jx)=pwaveH(i,j)
                                endif
                                keeploop = 1
                            endif
                        endif
                    endif
                end do
            endif
        end do
      end do
      if(keeploop == 0)exit
  end do

  pwaveH = pwaveH * ks

  return

end subroutine airymodel

subroutine wavtransport(pyits,pydepth,pyhent,pytransX,pytransY,pydz,pydist,pynrow,pyncol)

  implicit none

  integer :: pynrow,pyncol

  integer,intent(in) :: pyits
  real(kind=8),dimension(pynrow,pyncol),intent(in) :: pydepth
  real(kind=8),dimension(pynrow,pyncol),intent(in) :: pyhent
  real(kind=8),dimension(pynrow,pyncol),intent(in) :: pytransX
  real(kind=8),dimension(pynrow,pyncol),intent(in) :: pytransY

  real(kind=8),dimension(pynrow,pyncol),intent(out) :: pydz
  real(kind=8),dimension(pynrow,pyncol),intent(out) :: pydist

  integer :: i, j, k, loop, it, steps

  real(kind=8),dimension(pynrow,pyncol) :: ent,ndepth

  pydz = 0.
  steps = 20
  ndepth = pydepth+pyhent

  do k = 1, steps
    ent = pyhent/float(steps)
    loop = 0
    it = 0
    do while(loop==0 .and. it<pyits)
      loop = 1
      it = it+1
      do j = 2, pyncol-1
        do i = 2, pynrow-1
          if(ent(i,j)>0.01)then
            loop = 0
            ! Below critical shear stress for entrainment deposit everything
            if(pyhent(i,j)==0.)then
              pydz(i,j) = pydz(i,j)+ent(i,j)
              ndepth(i,j) = ndepth(i,j)-ent(i,j)
            else
              ! Along the X-axis

              ! Moving towards East
              if(pytransX(i,j)>0)then
                ! Inland deposit inside cell
                if(ndepth(i+1,j)<=0)then
                  pydz(i,j) = pydz(i,j)+pytransX(i,j)*ent(i,j)
                  ndepth(i,j) = ndepth(i,j)-pytransX(i,j)*ent(i,j)
                ! Transfert entrained sediment to neighbouring cell
                else
                  ! In case the directions are following the same trend
                  if(pytransX(i+1,j)>=0)then
                    ent(i+1,j) = ent(i+1,j)+pytransX(i,j)*ent(i,j)
                  ! In case the directions are facing each others
                  else
                    pydz(i,j) = pydz(i,j)+0.5*pytransX(i,j)*ent(i,j)
                    ndepth(i,j) = ndepth(i,j)-0.5*pytransX(i,j)*ent(i,j)
                    pydz(i+1,j) = pydz(i+1,j)+0.5*pytransX(i,j)*ent(i,j)
                    ndepth(i+1,j) = ndepth(i+1,j)-0.5*pytransX(i,j)*ent(i,j)
                  endif
                endif
              ! Moving towards West
              elseif(pytransX(i,j)<0)then
                ! Inland deposit inside cell
                if(ndepth(i-1,j)<=0)then
                  pydz(i,j) = pydz(i,j)-pytransX(i,j)*ent(i,j)
                  ndepth(i,j) = ndepth(i,j)+pytransX(i,j)*ent(i,j)
                ! Transfert entrained sediment to neighbouring cell
                else
                  ! In case the directions are following the same trend
                  if(pytransX(i-1,j)<=0)then
                    ent(i-1,j) = ent(i-1,j)-pytransX(i,j)*ent(i,j)
                  ! In case the directions are facing each others
                  else
                    pydz(i,j) = pydz(i,j)-0.5*pytransX(i,j)*ent(i,j)
                    ndepth(i,j) = ndepth(i,j)+0.5*pytransX(i,j)*ent(i,j)
                    pydz(i-1,j) = pydz(i-1,j)-0.5*pytransX(i,j)*ent(i,j)
                    ndepth(i-1,j) = ndepth(i-1,j)+0.5*pytransX(i,j)*ent(i,j)
                  endif
                endif
              endif

              ! Along the Y-axis

              ! Moving towards North
              if(pytransY(i,j)>0)then
                ! Inland deposit inside cell
                if(ndepth(i,j+1)<=0)then
                  pydz(i,j) = pydz(i,j)+pytransY(i,j)*ent(i,j)
                  ndepth(i,j) = ndepth(i,j)-pytransY(i,j)*ent(i,j)
                ! Transfert entrained sediment to neighbouring cell
                else
                  ! In case the directions are following the same trend
                  if(pytransY(i,j+1)>=0)then
                    ent(i,j+1) = ent(i,j+1)+pytransY(i,j)*ent(i,j)
                  ! In case the directions are facing each others
                  else
                    pydz(i,j) = pydz(i,j)+0.5*pytransY(i,j)*ent(i,j)
                    ndepth(i,j) = ndepth(i,j)-0.5*pytransY(i,j)*ent(i,j)
                    pydz(i,j+1) = pydz(i,j+1)+0.5*pytransY(i,j)*ent(i,j)
                    ndepth(i,j+1) = ndepth(i,j+1)-0.5*pytransY(i,j)*ent(i,j)
                  endif
                endif
              ! Moving towards South
              elseif(pytransY(i,j)<0)then
                ! Inland deposit inside cell
                if(ndepth(i,j-1)<=0)then
                  pydz(i,j) = pydz(i,j)-pytransY(i,j)*ent(i,j)
                  ndepth(i,j) = ndepth(i,j)+pytransY(i,j)*ent(i,j)
                ! Transfert entrained sediment to neighbouring cell
                else
                  ! In case the directions are following the same trend
                  if(pytransY(i,j-1)<=0)then
                    ent(i,j-1) = ent(i,j-1)-pytransY(i,j)*ent(i,j)
                  ! In case the directions are facing each others
                  else
                    pydz(i,j) = pydz(i,j)-0.5*pytransY(i,j)*ent(i,j)
                    ndepth(i,j) = ndepth(i,j)+0.5*pytransY(i,j)*ent(i,j)
                    pydz(i,j-1) = pydz(i,j-1)-0.5*pytransY(i,j)*ent(i,j)
                    ndepth(i,j-1) = ndepth(i,j-1)+0.5*pytransY(i,j)*ent(i,j)
                  endif
                endif

              endif
            endif
            ent(i,j) = 0.
          else
            pydz(i,j) = pydz(i,j)+ent(i,j)
            ndepth(i,j) = ndepth(i,j)+ent(i,j)
            ent(i,j) = 0.
          endif
        enddo
      enddo
    enddo
    if(it>=pyits)then
      pydz = pydz+ent
      ndepth = ndepth-ent
    endif
  enddo

  ! Find reworked sediment above water level
  pydist = 0.
  do j = 1, pyncol
    do i = 1, pynrow
      if(pydz(i,j)>pydepth(i,j)+pyhent(i,j).and.pydepth(i,j)+pyhent(i,j)>0.)then
        pydist(i,j) = pydz(i,j)-pydepth(i,j)-pyhent(i,j)
        pydz(i,j) = pydepth(i,j)+pyhent(i,j)
      endif
    enddo
  enddo

  return

end subroutine wavtransport

subroutine wavdiffusion(oelevpy,dzpy,coeffpy,maxthpy,tsteppy,nsteppy,depopy,numrowpy,numcolpy)

  implicit none
  
  integer :: numrowpy,numcolpy

  integer,intent(in) :: nsteppy
  real(kind=8),intent(in) :: maxthpy
  real(kind=8),intent(in) :: tsteppy
  real(kind=8),dimension(numrowpy,numcolpy),intent(in) :: oelevpy
  real(kind=8),dimension(numrowpy,numcolpy),intent(in) :: dzpy
  real(kind=8),intent(in) :: coeffpy


  real(kind=8),dimension(numrowpy,numcolpy),intent(out) :: depopy

  integer :: i,j,i1,j1,k,it
  real(kind=8),dimension(numrowpy,numcolpy) :: diffmarine,elev

  integer,dimension(4)::is=(/1,0,-1,0/)
  integer,dimension(4)::js=(/0,1,0,-1/)
  real(kind=8) :: flx,mindt

  depopy = dzpy
  elev = oelevpy

  do it=1,nsteppy
    diffmarine = 0.
    mindt = tsteppy
    do j = 2, numcolpy-1
      do i = 2, numrowpy-1
        do k = 1,4
          i1 = i+is(k)
          j1 = j+js(k)
          flx = elev(i1,j1)-elev(i,j)
          if(depopy(i,j)>maxthpy .and. elev(i,j)>elev(i1,j1))then
            diffmarine(i,j) = diffmarine(i,j) + flx*coeffpy
          elseif(depopy(i1,j1)>maxthpy .and. elev(i,j)<elev(i1,j1))then
            diffmarine(i,j) = diffmarine(i,j) + flx*coeffpy
          endif
        enddo
        if(diffmarine(i,j)<0. .and. diffmarine(i,j)*tsteppy<-depopy(i,j))then
          mindt = min(-depopy(i,j)/diffmarine(i,j),mindt)
        endif
      enddo
    enddo
    depopy = depopy + diffmarine*mindt
    elev = elev + diffmarine*mindt
  enddo

  return

end subroutine wavdiffusion
