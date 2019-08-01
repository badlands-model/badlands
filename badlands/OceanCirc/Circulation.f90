! =====================================================================================
! BADLANDS (BAsin anD LANdscape DynamicS)
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the icircmlied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Tecircmle
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  Circulation.f90
!
!    Description:  Ocean circulation model.
!
!        Version:  1.0
!        Created:  01/10/15 14:07:38
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module ocean_circ

  use parallel
  use bilinear
  use topology
  use parameters
  use external_forces

  implicit none


  integer::circn,circm,ocirc_steps

  integer,parameter::gap=2

  ! Conversion from m/s to knots
  real,parameter::knots=1.9438445
  real,parameter::grav=981.

  integer::xlen,ylen,xclen,yclen
  real,dimension(:),allocatable::x1_array,y1_array,oceanX,oceanY
  real,dimension(:),allocatable::xc_array,yc_array
  real,dimension(:,:),allocatable::o1,o2

  real::g1,g2,g3,g4,g5,g6,g7,alin
  real::oc_beta,oc_gamma,circdx,ah,atf,rho,ome,ocf,tiner,niner

  ! Extension values for circulation grid cocircmutational values
  integer::nn,mm,n,m,nnn,mmm,ne,nnnn,mmmm

  ! Circulation grid cocircmutational arrays
  integer,dimension(:),allocatable::ina,inb,jnab,imab,jma,jmb
  integer,dimension(:),allocatable::in1,in2,jn,im,jm1,jm2
  integer,dimension(:),allocatable::inn1,inn2,jnn,imm,jmm1,jmm2
  integer,dimension(:),allocatable::ie,je,ije

  integer::pt1,pt2,pt3,pt4,pt5,pt6,pt7

  real,dimension(:,:),allocatable::circ_a,circ_aa,circ_p
  integer,dimension(:,:),allocatable::p1,p2,p3,p4,p5,p6,p7

contains

  ! ============================================================================
  integer function ip(L)
    integer::L
    ip=int(L/2)*2
    return
  end function ip
  ! =====================================================================================
  subroutine circulation_compute(sgp)

    integer::kt,sgp
    real::akate,afilt,bfilt,tpp

    if(circON==0) return

    if(pet_id/=0) goto 30
    tpp=20.

    ! First cocircmute the ocean circulation
    do kt=1,ocirc_steps
       akate=(kt/real(ocirc_steps))*100.
       alin=1.
       alin=real(kt)/niner
       if(kt>=niner) alin=1.
       if(akate>=tpp)then
          if(pet_id==0) write(*,13) int(akate)
          tpp=tpp+20.
       endif
       if(kt==ocirc_steps-1) call record_circulation_velocity
       call update_circulation_bounds(kt)
       afilt=real(kt)/real(oceanfilter)
       bfilt=real(kt/oceanfilter)
       if(afilt==oceanfilter) call circulation_filter
       call sweepx_fct
       call update_circulation_bounds(kt)
       call sweepy_fct
    enddo

13  format('+',' currents circulation completion: ',i5,' %')

30  continue

    call mpi_barrier(badlands_world,rc)
    call get_circulation_velocity(sgp)

    return

  end subroutine circulation_compute
  ! =====================================================================================
  subroutine get_circulation_velocity(sgp)

    integer::sgp,i,j,p
    real(kind=8),dimension(nbox*nboy)::velcircV,velcircU

    if(.not.allocated(rcircU))then
      allocate(rcircU(nbox,nboy))
      allocate(rcircV(nbox,nboy))
    else
      p=0
      do i=1,nbox
        do j=1,nboy
          p=p+1
          velcircV(p)=real(rcircV(i,j),8)
          velcircU(p)=real(rcircU(i,j),8)
        enddo
      enddo
    endif

    if(pet_id/=0)then
      rcircU=-1000.
      rcircV=-1000.
    endif

    call mpi_bcast(velcircV,nbox*nboy,mpi_double_precision,0,badlands_world,rc)
    call mpi_bcast(velcircU,nbox*nboy,mpi_double_precision,0,badlands_world,rc)

    p=0
    do i=1,nbox
      do j=1,nboy
        p=p+1
        circU(i,j,sgp)=real(velcircU(p)/100.)
        circV(i,j,sgp)=real(velcircV(p)/100.)
      enddo
    enddo

  end subroutine get_circulation_velocity
  ! =====================================================================================
  subroutine record_circulation_velocity

    integer::k,i,i1,i2,j1,j
    real::x1,y1
    real,dimension(circm,circn)::Umean,Vmean

    if(.not.allocated(rcircU))then
      allocate(rcircU(nbox,nboy))
      allocate(rcircV(nbox,nboy))
    endif
    Umean=0.
    Vmean=0.

    do i=1,circm
      do j=1,circn
        circ_aa(i,j)= 0.0
      enddo
    enddo

    do k=1,nn
      i1=in1(k)
      i2=in2(k)
      j=jn(k)
      if(i1==ip(i1)) i1=i1+1
      if(i2==ip(i2)) i2=i2-1
      do i=i1,i2,2
        if(i>=1.and.i<=circm) circ_aa(i,j)=circ_a(i,j)
      enddo
    enddo

    ! Get velocities and directions
    do k=1,nn
      i1=in1(k)
      i2=in2(k)
      j=jn(k)
      if(i1/=ip(i1)) i1=i1-1
      if(i2/=ip(i2)) i2=i2+1
      do i=i1,i2,2
        if(i>=1.and.i<=circm)then
          Vmean(i,j)=circ_a(i,j)
        endif
        if(i>1.and.i<circm)then
          Umean(i,j)=0.25*(circ_a(i+1,j+1)+circ_a(i+1,j-1)+circ_a(i-1,j-1)+circ_a(i-1,j+1))
        endif
      enddo
    enddo

    ! Sort computed values points
    i1=0
    j1=0
    do i=1,circm
      if(i==ip(i)) i1=i1+1
      do j=1,circn
        if(i==ip(i).and.j==ip(j))then
          j1=j1+1
          o1(j1,i1)=Umean(i,j)
          o2(j1,i1)=Vmean(i,j)
        endif
      enddo
      j1=0
    enddo

    ! Interpolate velocities over the computational grid
    do i=1,circm
      do j=1,circn
        if((j-1)*ocean_dx>x1_array(1).and.(j-1)*ocean_dx<x1_array(xlen))then
          if((i-1)*ocean_dx>y1_array(1).and.(i-1)*ocean_dx<y1_array(ylen))then
            y1=(i-1)*ocean_dx
            x1=(j-1)*ocean_dx
            Umean(i,j)=interpolate_circulation_grid(xlen,x1_array(1:xlen),ylen,y1_array(1:ylen),o1,x1,y1)
            Vmean(i,j)=interpolate_circulation_grid(xlen,x1_array(1:xlen),ylen,y1_array(1:ylen),o2,x1,y1)
          endif
        endif
      enddo
    enddo

    ! Interpolate velocities over stratigraphic mesh
    do i=1,nbox
      do j=1,nboy
        if(oceanZ(i,j)<0.0.and.oceanZ(i,j)>-circbase)then
          rcircU(i,j)=Umean(i+gap,j+gap)
          rcircV(i,j)=Vmean(i+gap,j+gap)
        else
          rcircU(i,j)=0.0
          rcircV(i,j)=0.0
        endif
      enddo
    enddo

    rcircU(1,:)=rcircU(3,:)
    rcircU(2,:)=rcircU(3,:)
    rcircU(nbox,:)=rcircU(nbox-2,:)
    rcircU(nbox-1,:)=rcircU(nbox-2,:)
    rcircU(:,1)=rcircU(:,3)
    rcircU(:,2)=rcircU(:,3)
    rcircU(:,nboy)=rcircU(:,nboy-2)
    rcircU(:,nboy-1)=rcircU(:,nboy-2)

    rcircV(1,:)=rcircV(3,:)
    rcircV(2,:)=rcircV(3,:)
    rcircV(nbox,:)=rcircV(nbox-2,:)
    rcircV(nbox-1,:)=rcircV(nbox-2,:)
    rcircV(:,1)=rcircV(:,3)
    rcircV(:,2)=rcircV(:,3)
    rcircV(:,nboy)=rcircV(:,nboy-2)
    rcircV(:,nboy-1)=rcircV(:,nboy-2)

    do i=1,circm
      do j=1,circn
        if(i==ip(i).and.j/=ip(j)) goto 500
        circ_aa(i,j)=circ_a(i,j)
500     continue
      enddo
    enddo

    return

  end subroutine record_circulation_velocity
  ! =====================================================================================
  subroutine update_circulation_bounds(kt)

    integer::kt,km,kk,k1,k2,kn,k,halo

    real::amu1,zcrit,amu,amumu

    if(kt>=niner) goto 7
    halo=2

    do km=1,ne
       if(ie(km)<=halo+1)then
          if(kt==1) circ_a(ie(km),je(km))=(float(kt)/niner)*circ_a(ie(km),je(km))
          circ_aa(ie(km),je(km))=(float(kt)/niner)*circ_aa(ie(km),je(km))
       elseif(ie(km)>=circm-halo+1)then
          if(kt==1) circ_a(ie(km),je(km))=(float(kt)/niner)*circ_a(ie(km),je(km))
          circ_aa(ie(km),je(km))=(float(kt)/niner)*circ_aa(ie(km),je(km))
       endif
    enddo

    ! Solves the partially clacircmed explicit radiation condition
    if(kt==1) goto 7

    amu1=ocean_dt/atf
    zcrit=4427.0
    do km=1,ne
      ! Bottom margin
      if(ije(km)/=-1) goto 13
      if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)<circm.and.je(km)>1.and.je(km)<circn) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)+1,je(km)-1)+circ_a(ie(km)+1,je(km)+1)))/circdx
       if(ie(km)==ip(ie(km)).and.je(km)==ip(je(km)).and.je(km)>1.and.je(km)<circn.and.ie(km)>0) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km),je(km)-1)+circ_a(ie(km),je(km)+1)))/circdx
       if(ie(km)>0.and.je(km)>0.and.ie(km)+4<circm)then
          circ_aa(ie(km),je(km))=(circ_a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*circ_a(ie(km)+2,je(km))-circ_a(ie(km)+4,je(km))- &
               3.0*circ_a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*circdx/ocean_dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) circ_aa(ie(km),je(km))=0.0

       ! Upper margin
13     if(ije(km)/=1) goto 14
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.je(km)>1.and.je(km)<circn) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)-1,je(km)-1)+circ_a(ie(km)-1,je(km)+1)))/circdx
       if(ie(km)==ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>0.and.je(km)>1.and.je(km)<circn) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km),je(km)-1)+circ_a(ie(km),je(km)+1)))/circdx
       if(ie(km)<circn.and.ie(km)>0.and.je(km)>0.and.ie(km)-4>0)then
          circ_aa(ie(km),je(km))=(circ_a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*circ_a(ie(km)-2,je(km))-circ_a(ie(km)-4,je(km))- &
               3.0*circ_a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*circdx/ocean_dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) circ_aa(ie(km),je(km))=0.0

       ! Left margin
14     if(ije(km)/=-2) goto 15
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.ie(km)<circm.and.je(km)<circn) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)+1,je(km)+1)+circ_a(ie(km)-1,je(km)+1)))/circdx
       if(ie(km)/=ip(ie(km)).and.je(km)/=ip(je(km)).and.ie(km)>1.and.ie(km)<circm.and.je(km)>0) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)+1,je(km))+circ_a(ie(km)-1,je(km))))/circdx
       if(ie(km)>0.and.je(km)>0.and.je(km)+4<circn)then
          circ_aa(ie(km),je(km))=(circ_a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*circ_a(ie(km),je(km)+2)-circ_a(ie(km),je(km)+4)- &
               3.0*circ_a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*circdx/ocean_dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) circ_aa(ie(km),je(km))=0.0

       ! Right margin
15     if(ije(km)/=2) goto 11
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.ie(km)<circm.and.je(km)>1) &
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)+1,je(km)-1)+circ_a(ie(km)-1,je(km)-1)))/circdx
       if(ie(km)/=ip(ie(km)).and.je(km)/=ip(je(km)).and.ie(km)>1.and.ie(km)<circm.and.je(km)>0)&
            amu=ocean_dt*sqrt(grav*0.5*(circ_a(ie(km)+1,je(km))+circ_a(ie(km)-1,je(km))))/circdx
       if(ie(km)>0.and.je(km)>0.and.je(km)-4>0)then
          circ_aa(ie(km),je(km))=(circ_a(ie(km),je(km))*(1.0-0.25*amu1)+&
               amu*0.125*(4.0*circ_a(ie(km),je(km)-2)-circ_a(ie(km),je(km)-4)- &
               3.0*circ_a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*circdx/ocean_dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) circ_aa(ie(km),je(km))=0.0

11     continue
    enddo

7   do kk=1,nnn
       if(ina(kk)>0.and.ina(kk)+2<=circm.and.jnab(kk)>0)&
            circ_a(ina(kk),jnab(kk))=circ_a(ina(kk)+2,jnab(kk))
      if(inb(kk)>2.and.inb(kk)<circm.and.jnab(kk)>0)&
            circ_a(inb(kk),jnab(kk))=circ_a(inb(kk)-2,jnab(kk))
    enddo

    do kk=1,mmm
       if(imab(kk)>0.and.jma(kk)>0.and.jma(kk)+2<=circn)&
            circ_a(imab(kk),jma(kk))=circ_a(imab(kk),jma(kk)+2)
       if(imab(kk)>0.and.jmb(kk)>2)&
            circ_a(imab(kk),jmb(kk))=circ_a(imab(kk),jmb(kk)-2)
    enddo

    do k1=1,mmm
       do k2=1,nnn
          if(ina(k2)==imab(k1).and.jnab(k2)==jma(k1).and.ina(k2)>0.and.jnab(k2)>0.and.ina(k2)+2<=circm.and.jnab(k2)+2<=circn) &
               circ_a(ina(k2),jnab(k2))=(circ_a(ina(k2)+2,jnab(k2))+circ_a(ina(k2),jnab(k2)+2))/2.
          if(inb(k2)==imab(k1).and.jnab(k2)==jma(k1).and.inb(k2)>2.and.jnab(k2)>0.and.jnab(k2)+2<=circn) &
               circ_a(inb(k2),jnab(k2))=(circ_a(inb(k2)-2,jnab(k2))+circ_a(inb(k2),jnab(k2)+2))/ 2.
          if(ina(k2)==imab(k1).and.jnab(k2)==jmb(k1).and.ina(k2)>0.and.ina(k2)+2<=circm.and.jnab(k2)>2) &
               circ_a(ina(k2),jnab(k2))=(circ_a(ina(k2),jnab(k2)-2)+circ_a(ina(k2)+2,jnab(k2)))/2.
          if(inb(k2)==imab(k1).and.jnab(k2)==jmb(k1).and.inb(k2)>2.and.jnab(k2)>2) &
               circ_a(inb(k2),jnab(k2))=(circ_a(inb(k2)-2,jnab(k2))+circ_a(inb(k2),jnab(k2)-2))/2.
       enddo
    enddo

    do kn=1,nn
       do k1=1,ne
          if(in1(kn)==ie(k1).and.jn(kn)==je(k1)) goto 60
       enddo
       if(in1(kn)>0.and.jn(kn)>0.and.in1(kn)<=circm.and.jn(kn)<=circn)then
          circ_a(in1(kn),jn(kn))=0.0
          circ_aa(in1(kn),jn(kn))=0.0
       endif
60     do k1=1,ne
          if(in2(kn)==ie(k1).and.jn(kn)==je(k1)) goto 30
       enddo
       if(in2(kn)>0.and.jn(kn)>0.and.in2(kn)<=circm.and.jn(kn)<=circn)then
          circ_a(in2(kn),jn(kn))=0.0
          circ_aa(in2(kn),jn(kn))=0.0
       endif
30     continue
    enddo

    do km=1,mm
       do k1=1,ne
          if(im(km)==ie(k1).and.jm1(km)==je(k1)) goto 90
       enddo
       if(im(km)>0.and.jm1(km)>0.and.im(km)<=circm.and.jm1(km)<=circn)then
          circ_a(im(km),jm1(km))=0.0
          circ_aa(im(km),jm1(km))=0.0
       endif
90     do k1=1,ne
          if(im(km)==ie(k1).and.jm2(km)==je(k1)) goto 40
       enddo
       if(im(km)>0.and.jm2(km)>0.and.im(km)<=circm.and.jm2(km)<=circn)then
          circ_a(im(km),jm2(km))=0.0
          circ_aa(im(km),jm2(km))=0.0
       endif
40     continue
    enddo

    do k=1,nnnn
       if(inn1(k)/=ip(inn1(k)).and.inn1(k)>0.and.jnn(k)>0.and.inn1(k)+2<=circm.and.jnn(k)<=circn) &
            circ_a(inn1(k),jnn(k))=circ_a(inn1(k)+2,jnn(k))
       if(inn2(k)/=ip(inn2(k)).and.inn2(k)>2.and.jnn(k)>0.and.inn2(k)<=circm.and.jnn(k)<=circn) &
            circ_a(inn2(k),jnn(k))=circ_a(inn2(k)-2,jnn(k))
    enddo

    do k=1,mmmm
      if(jmm1(k)==ip(jmm1(k)).and.imm(k)>0.and.jmm1(k)>0.and.imm(k)<=circm.and.jmm1(k)+2<=circn) &
        circ_a(imm(k),jmm1(k))=circ_a(imm(k),jmm1(k)+2)
      if(jmm2(k)==ip(jmm2(k)).and.imm(k)>0.and.jmm2(k)>2.and.imm(k)<=circm.and.jmm2(k)<=circn) &
        circ_a(imm(k),jmm2(k))=circ_a(imm(k),jmm2(k)-2)
    enddo

    return

  end subroutine update_circulation_bounds
  ! =====================================================================================
  subroutine sweepx_fct

    integer::kn,i,j,i1,i2,ide,idid

    real::X

    real::U1(0:circm),U2(0:circm),W1(0:circm),W2(0:circm)

    U1=0.
    U2=0.
    W1=0.
    W2=0.

    do kn=1,nn
       j=jn(kn)
       i1=in1(kn)
       i2=in2(kn)-2
       if(ip(i1)/=i1.and.ip(i2)==i2) ide=1
       if(ip(i1)==i1.and.ip(i2)==i2) ide=2
       if(ip(i1)==i1.and.ip(i2)/=i2) ide=3
       if(ip(i1)/=i1.and.ip(i2)/=i2) ide=4

       if(ide==1.or.ide==3) i1=i1+1

       lp: do  i=i1,i2,2

          if(i > i2) exit lp

          idid=0
          if(ide==2.and.i==i1) idid=1
          if(i>2.and.i<circm-2)then
             call gegeX(i,j,ide)
          else
             goto 40
          endif
          if(ide==2.and.i==i1) goto 50
          if(ide==4.and.i==i1) goto 60
          if(ide==3.or.ide==4) goto 45

          if(i==i1) U2(i-2)=circ_aa(i-1,j)
          if(i==i1) W2(i-2)=0.0
          X=g1+g2*W2(i-2)-g3*g5
          U1(i)=(g4+g2*U2(i-2)-g3*g7)/X
          W1(i)=g3*g6/X
          U2(i)=(g7*(g1+g2*W2(i-2))-g5*(g4+g2*U2(i-2)))/X
          W2(i)=g6*(g1+g2*W2(i-2))/X
          goto 40

45        if(i==i1) U1(i-2)=circ_aa(i-1,j)
          if(i==i1) W1(i-2)=0.0
          X=g1*(1-g5*W1(i-2))+g2*g6
          U1(i)=(g2*(g7-g5*U1(i-2))+g4*(1-g5*W1(i-2)))/X
          W1(i)=g3*(1-g5*W1(i-2))/X
          U2(i)=(g1*(g7-g5*U1(i-2))-g6*g4)/X
          W2(i)=g3*g6/X
          goto 40

50        U2(i)=g7-g5*circ_aa(i,j)
          W2(i)=g6
          U1(i)=0.0
          W1(i)=0.0
          goto 40

60        X= g1
          U1(i)=(g4+g2*circ_aa(i,j))/X
          W1(i)=g3/X
40        continue
       enddo lp

       !  Computes UU and EE values
       lp2: do i=i2,i1,-2
          if(i < i1) exit lp2
          if(i2>circm-2) goto 70
          if(ide==4.and.i==i1) goto 80
          if(ide==2.and.i==i1) goto 90
          if(ide==3.or.ide==4) goto 85

          circ_aa(i,j)=U1(i)+W1(i)*circ_aa(i+2,j)
          circ_aa(i+1,j)=U2(i)-W2(i)*circ_aa(i+2,j)
          goto 70
85        circ_aa(i+1,j)=U1(i)-W1(i)*circ_aa(i+2,j)
          circ_aa(i,j)=U2(i)+W2(i)*circ_aa(i+2,j)
          goto 70
80        circ_aa(i+1,j)=U1(i)-W1(i)*circ_aa(i+2,j)
          goto 70
90        circ_aa(i+1,j)=U2(i)-W2(i)*circ_aa(i+2,j)
70        continue

       enddo lp2
    enddo
    call calculateV

    do i=1,circm
       do j=1,circn
          if(i==ip(i).and.j/=ip(j)) goto 150
          circ_a(i,j)=circ_aa(i,j)
150       continue
       enddo
    enddo

    return

  end subroutine sweepx_fct
  ! =====================================================================================
  subroutine sweepy_fct

    integer::kn,i,j,j1,j2,ide,idid

    real::X

    real::U1(0:circn),U2(0:circn),W1(0:circn),W2(0:circn)

    U1=0.0
    U2=0.0
    W1=0.0
    W2=0.0

    do kn=1,mm

       i=im(kn)
       j1=jm1(kn)
       j2=jm2(kn)-2

       if(ip(j1)==j1.and.ip(j2)/=j2) ide=1
       if(ip(j1)/=j1.and.ip(j2)/=j2) ide=2
       if(ip(j1)/=j1.and.ip(j2)==j2) ide=3
       if(ip(j1)==j1.and.ip(j2)==j2) ide=4

       if(ide==1.or.ide==3) j1=j1+1

       do j=j1,j2,2

          idid=0
          if(ide==2.and.j==j1) idid=1

          if(i>2.and.i<circm-2)then
             call gegeY(i,j,ide)
          else
             goto 40
          endif
          if(ide==2.and.j==j1) goto 50
          if(ide==4.and.j==j1) goto 60
          if(ide==3.or.ide==4) goto 45

          if(j==j1) U2(j-2)=circ_aa(i,j-1)
          if(j==j1) W2(j-2)=0.

          X=g1+g2*W2(j-2)-g3*g5
          U1(j)=(g4+g2*U2(j-2)-g3*g7)/X
          W1(j)=g3*g6/X
          U2(j)=(g7*(g1+g2*W2(j-2))-g5*(g4+g2*U2(j-2)))/X
          W2(j)=g6*(g1+g2*W2(j-2))/X
          goto 40

45        if(j==j1) U1(j-2)=circ_aa(i,j-1)
          if(j==j1) W1(j-2)=0.
          X=g1*(1-g5*W1(j-2))+g2*g6
          U1(j)=(g2*(g7-g5*U1(j-2))+g4*(1-g5*W1(j-2)))/X
          W1(j)=g3*(1.0-g5*W1(j-2))/X
          U2(j)=(g1*(g7-g5*U1(j-2))-g6*g4)/X
          W2(j)=g3*g6/X
          goto 40

50        U2(j)=g7-g5*circ_aa(i,j)
          W2(j)=g6
          U1(j)=0.
          W1(j)=0.
          goto 40

60        X=g1
          U1(j)=(g4+g2*circ_aa(i,j))/X
          W1(j)=g3/X
40        continue
       enddo

       ! Cocircmutes the VV and ee values
       do j=j2,j1,-2

          if((ide==4).and.(j==j1)) goto 80
          if((ide==2).and.(j==j1)) goto 90
          if((ide==3).or. (ide==4)) goto 85

          circ_aa(i,j)=U1(j)+W1(j)*circ_aa(i,j+2)
          circ_aa(i,j+1)=U2(j)-W2(j)*circ_aa(i,j+2)
          goto 70
85        circ_aa(i,j+1)=U1(j)-W1(j)*circ_aa(i,j+2)
          circ_aa(i,j)=U2(j)+W2(j)*circ_aa(i,j+2)
          goto 70
80        circ_aa(i,j+1)=U1(j)-W1(j)*circ_aa(i,j+2)
          goto 70
90        circ_aa(i,j+1)=U2(j)-W2(j)*circ_aa(i,j+2)
70        continue
       enddo
    enddo
    call calculateU

    do i=1,circm
       do j=1,circn
          if(i==ip(i).and.j/=ip(j)) goto 150
          circ_a(i,j)=circ_aa(i,j)
150       continue
       enddo
    enddo

    return

  end subroutine sweepy_fct
  ! =====================================================================================
  subroutine circulation_filter

    integer::i,j,k,i1,i2,j1,j2
    real::S

    S=0.1

    do k=1,nn
       i1=in1(k)
       i2=in2(k)
       if(ip(i1)/=i1) i1=i1+1
       if(ip(i2)/=i2) i2=i2-1
       j=jn(k)
       do i=i1,i2,2
          if(i>2.and.i<=circm-2)then
             circ_aa(i,j)=circ_a(i,j)+(S/2.0)*(1.0-S)*(circ_a(i+2,j)+circ_a(i-2,j)+circ_a(i,j-2)+circ_a(i,j+2) &
                -4.0*circ_a(i,j))+(S*S/4.0)*(circ_a(i-2,j+2)+circ_a(i+2,j+2)+circ_a(i+2,j-2) &
                 +circ_a(i-2,j-2)-4.0*circ_a(i,j))
          endif
       enddo
    enddo

    do k=1,mm
       j1=jm1(k)
       j2=jm2(k)
       if(ip(j1)==j1) j1=j1+1
       if(ip(j2)==j2) j2=j2-1
       i=im(k)
       do j=j1,j2,2
          if(i>2.and.i<=circm-2)then
             circ_aa(i,j)=circ_a(i,j)+(S/2.0)*(1.0-S)*(circ_a(i+2,j)+circ_a(i-2,j)+circ_a(i,j-2)+circ_a(i,j+2) &
                -4.0*circ_a(i,j))+(S*S/4.0)*(circ_a(i-2,j+2)+circ_a(i+2,j+2)+circ_a(i+2,j-2) &
                 +circ_a(i-2,j-2)-4.0*circ_a(i,j))
          endif
       enddo
    enddo

    do k=1,nn
       i1=in1(k)+2
       i2=in2(k)-2
       if(ip(i1)==i1) i1=i1+1
       if(ip(i2)==i2) i2=i2-1
       j=jn(k)
       do i=i1,i2,2
          if(i>2.and.i<=circm-2)then
             circ_aa(i,j)=circ_a(i,j)+(S/2.0)*(1.0-S)*(circ_a(i+2,j)+circ_a(i-2,j)+circ_a(i,j-2)+circ_a(i,j+2) &
                -4.0*circ_a(i,j))+(S*S/4.0)*(circ_a(i-2,j+2)+circ_a(i+2,j+2)+circ_a(i+2,j-2) &
                 +circ_a(i-2,j-2)-4.0*circ_a(i,j))
          endif
       enddo
    enddo
    do i=1,circm
       do j=1,circn
          if(i==ip(i).and.j/=ip(j)) goto 70
          circ_a(i,j)=circ_aa(i,j)
70        continue
       enddo
    enddo

    return

  end subroutine circulation_filter
  ! =====================================================================================
  subroutine gegeY(i,j,ide)

    integer::i,j,ide,jaUXi

    real::S,b,dtx,ViXm,VimOD,Ca,uhm,vhm
    real::TSX,TSY,D,DDD,Um,f1,f2,f3,f4,f5,f6,f7,f8
    real::f9,f10,f11,f12,f13,f14,f15,f16,f17

    S=0.
    b=oc_beta
    dtx=ocean_dt/circdx

    if(ide==3.or.ide==4) then
       jaUXi=j
       j=j+1
    endif

    ViXm=0.25*(circ_p(i+1,j-1)+circ_p(i+1,j+1)+circ_p(i-1,j+1)+circ_p(i-1,j-1))

    VimOD=sqrt(circ_p(i,j)*circ_p(i,j)+ViXm*ViXm)
    VimOD=VimOD*51.44

    Ca=0.0012875
    if(VimOD>750.) Ca=0.0008+6.5e-07*VimOD

    TSX=0.001*Ca*51.44*ViXm*VimOD
    TSY=0.001*Ca*51.44*circ_p(i,j)*VimOD
    TSX=TSX*alin
    TSY=TSY*alin
    D=0.5*(circ_a(i-1,j)+circ_a(i+1,j)+circ_a(i,j+1)+circ_a(i,j-1))
    if(D<=0.) return

    DDD=0.5*(circ_a(i-1,j)+circ_a(i+1,j))

    Um=0.25*(circ_a(i-1,j-1)+circ_a(i+1,j-1)+circ_a(i-1,j+1)+circ_a(i+1,j+1))

    f1=0.5*S*circ_a(i,j)/D

    f2=0.125*f1*dtx*(circ_a(i+1,j+2)-circ_a(i+1,j-2)+circ_a(i-1,j+2) &
       -circ_a(i-1,j-2))

    f3=0.25*S*dtx*Um*(circ_a(i+1,j)-circ_a(i-1,j))/D

    f4=0.0625*S*dtx*Um*(circ_a(i+2,j+1)-circ_a(i-2,j+1)+&
         circ_a(i+2,j-1)-circ_a(i-2,j-1))/D

    f5=0.125*S*dtx*(circ_a(i+1,j+1)-circ_a(i-1,j+1)+&
         circ_a(i+1,j-1)-circ_a(i-1,j-1))

    f6=0.125*S*dtx*Um*(circ_a(i+2,j)-circ_a(i-2,j))
    f7=0.5*ocean_dt*ocf*Um
    f8=0.25*dtx*(grav-S*(circ_a(i,j)**2)/D)
    f9=0.5*ocean_dt*TSY/(rho*D)

    uhm=Um
    vhm=circ_a(i,j)

    f10=0.5*ocean_dt*oceanfric*sqrt(uhm*uhm+vhm*vhm)/D
    f11=0.25*dtx*oc_gamma*(circ_p(i,j+1)-circ_p(i,j-1))/(rho/1000.)
    f11=f11*alin
    f12=0.125*dtx*ah*(circ_a(i+2,j)-2*circ_a(i,j)+circ_a(i-2,j))/circdx
    f13=0.125*dtx*ah*(circ_a(i,j+2)-2*circ_a(i,j)+circ_a(i,j-2))/circdx

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j-1
    endif

    f14=0.0625*dtx*(circ_a(i+1,j+1)+circ_a(i-1,j+1))*(circ_a(i+1,j+2) &
       -circ_a(i-1,j+2)+circ_a(i+1,j)-circ_a(i-1,j)+circ_a(i+2,j+1)-&
         circ_a(i-2,j+1))

    f15=0.0625*dtx*(circ_a(i+1,j+2)-circ_a(i+1,j)+circ_a(i-1,j+2) &
       -circ_a(i-1,j)+circ_a(i,j+3)-circ_a(i,j-1))

    f16=0.0625 *dtx*(4*circ_a(i,j+1)+circ_a(i+1,j+2)+circ_a(i+1,j)+&
         circ_a(i-1,j)+circ_a(i-1,j+2))*(circ_a(i+1,j+1)-circ_a(i-1,j+1))

    f17=0.0625*dtx*(4*circ_a(i,j+1)+circ_a(i+1,j+2)+circ_a(i+1,j) &
        +circ_a(i-1,j)+circ_a(i-1,j+2))

    g1=1.0-b*(f2+f3+f4+f5-f10)
    g2=f1+b*f8
    g3=b*f8-f1

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j+1
    endif

    g4=circ_a(i,j)*(1.0+(1.0-b)*(f2+f3+f4+f5-f10))+&
         circ_a(i,j-1)*((1.0-b)*f8-f1)-circ_a(i,j+1)*((1.0-b)*f8+&
         f1)-f6-f7+f9-f11+f12+f13
    g5=b*(f15-f17)
    g6=b*(f15+f17)

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j-1
    endif

    g7=circ_a(i,j+1)+circ_a(i,j)*(1.0-b)*(f17-f15)-circ_a(i,j+2)*&
         (1.0-b)*(f15+f17)-f14-f16

    if(ide==3.or.ide==4) j=jaUXi

    return

  end subroutine gegeY
  ! =====================================================================================
  subroutine gegeX(i,j,ide)

    integer::i,j,ide,iaUXi

    real::S,b,dtx,ViYm,VimOD,Ca,uhm,vhm
    real::TSX,TSY,D,DDD,Vm,f1,f2,f3,f4,f5,f6,f7,f8
    real::f9,f10,f11,f112,f113,f12,f13,f14,f15

    S=0.
    b=oc_beta
    dtx=ocean_dt/circdx

    if(ide==3.or.ide==4)then
       iaUXi=i
       i=i+1
    endif

    ViYm=0.25*(circ_p(i+1,j-1)+circ_p(i+1,j+1)+circ_p(i-1,j+1)+circ_p(i-1,j-1))
    VimOD=sqrt(circ_p(i,j)*circ_p(i,j)+ViYm*ViYm)
    VimOD=VimOD*51.44

    Ca=0.0012875
    if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD

    TSX=0.001*Ca*51.44*circ_p(i,j)*VimOD
    TSY=0.001*Ca*51.44*ViYm*VimOD
    TSX=TSX*alin
    TSY=TSY*alin

    D=0.5*(circ_a(i,j-1)+circ_a(i,j+1)+circ_a(i+1,j)+circ_a(i-1,j))
    if(D <= 0.) return

    DDD=0.5*(circ_a(i,j-1)+circ_a(i,j+1))
    Vm=0.25*(circ_a(i-1,j-1)+circ_a(i+1,j-1)+circ_a(i-1,j+1)+circ_a(i+1,j+1))

    f1=0.5*S*circ_a(i,j)/D

    f2=0.125*f1*dtx*(circ_a(i+2,j-1)-circ_a(i-2,j-1)+circ_a(i+2,j+1) &
       -circ_a(i-2,j+1))
    f3=0.25*S*dtx*Vm*(circ_a(i,j+1)-circ_a(i,j-1))/D

    f4=0.0625*S*dtx*Vm*(circ_a(i+1,j+2)-circ_a(i+1,j-2)+&
         circ_a(i-1,j+2)-circ_a(i-1,j-2))/D
    f5=0.125*S*dtx*Vm*(circ_a(i,j+2)-circ_a(i,j-2))

    f6=0.125*S*dtx *(circ_a(i+1,j+1)-circ_a(i+1,j-1)+circ_a(i-1,j+1) &
       -circ_a(i-1,j-1))

    f7=0.5*ocean_dt*ocf*Vm
    f8=0.25*dtx*(grav-S*((circ_a(i,j))**2)/D)
    f9=0.5*ocean_dt*TSX/(rho*D)

    uhm=circ_a(i,j)
    vhm=Vm

    f10=0.5*ocean_dt*oceanfric*sqrt(uhm*uhm+vhm*vhm)/D
    f11=0.25*dtx*oc_gamma*(circ_p(i+1,j)-circ_p(i-1,j))/(rho/1000.)
    f11=f11*alin
    f112=0.125*dtx*ah*(circ_a(i+2,j)-2.0*circ_a(i,j)+circ_a(i-2,j))/circdx
    f113=0.125*dtx*ah*(circ_a(i,j+2)-2.0*circ_a(i,j)+circ_a(i,j-2))/circdx

    if(ide==3.or.ide==4)then
       i=iaUXi
       i=i-1
    endif

    f12=0.0625*dtx*(circ_a(i+2,j-1)-circ_a(i,j-1)+circ_a(i+2,j+1) &
       -circ_a(i,j+1)+circ_a(i+3,j)-circ_a(i-1,j))

    f13=0.0625*dtx*(4*circ_a(i+1,j)+circ_a(i,j+1)+circ_a(i,j-1) &
        +circ_a(i+2,j+1)+circ_a(i+2,j-1))

    f14=0.0625*dtx*(circ_a(i+1,j+1)+circ_a(i+1,j-1))*(circ_a(i+2,j+1) &
       -circ_a(i+2,j-1)+circ_a(i,j+1)-circ_a(i,j-1)+circ_a(i+1,j+2)-&
         circ_a(i+1,j-2))

    f15=0.0625*dtx*(4*circ_a(i+1,j)+circ_a(i+2,j+1)+circ_a(i+2,j-1) &
        +circ_a(i,j+1)+circ_a(i,j-1))*(circ_a(i+1,j+1)-circ_a(i+1,j-1))

    g1=1-b*(f2+f3+f4+f6-f10)
    g2=f1+b*f8
    g3=b*f8-f1

    if(ide==3.or.ide==4) then
       i=iaUXi
       i=i+1
    endif

    g4=circ_a(i,j)*(1.0+(1.0-b)*(f2+f3+f4+f6-f10))+&
         circ_a(i-1,j)*((1.0-b)*f8-f1)-circ_a(i+1,j)*((1.0-b)*f8+&
         f1)-f5+f7+f9-f11+f112+f113
    g5=b*(f12-f13)
    g6=b*(f12+f13)

    if(ide==3.or.ide==4) then
       i=iaUXi
       i=i-1
    endif

    g7=circ_a(i+1,j)+circ_a(i,j)*(1.0-b)*(f13-f12)-circ_a(i+2,j)*&
         (1.0-b)*(f12+f13)-f14-f15

    if(ide==3.or.ide==4) i=iaUXi

    return

  end subroutine gegeX
  ! =====================================================================================
  subroutine calculateV

    integer::j,i,kn,iSen1,iSen2,j1,j2,ideS1,ini,ifi

    real::S,b,dtx,Um,D,DDD,f1,f2,f3,f4,f5,f6,f7,f8,f9

    real::ViXm,VimOD,Ca,X,TSX,TSY,uhm,vhm

    real::U1(0:circn),W1(0:circn)

    U1=0.
    W1=0.

    dtx=ocean_dt/circdx
    S=0.
    b=oc_beta

    do kn=1,mm
       iSen1=0
       iSen2=0
       i=im(kn)
       if(i<=2.or.i>circm-2) goto 10
       j1= jm1(kn)+2
       j2= jm2(kn)-2
       ideS1=jm2(kn)-jm1(kn)
       if(ideS1==2) goto 10

       ini=j1
       ifi=j2
       if(ip(j1)/=j1) goto 20

       j1=j1-1
       ini=j1
       iSen1=1

20     if(ip(j2)/=j2) goto 30
       j2=j2+1
       ifi=j2
       iSen2=1

30     j=ini
35     if(iSen1==0.and.iSen2==1.and.j==ini) j=ifi

       Um=0.25*(b*(circ_aa(i+1,j+1)+circ_aa(i+1,j-1)+circ_aa(i-1,j-1)+&
            circ_aa(i-1,j+1))+(1.0-b)*(circ_a(i+1,j+1)+circ_a(i+1,j-1)+&
            circ_a(i-1,j-1)+circ_a(i-1,j+1)))
       D=0.5*(circ_a(i+1,j)+circ_a(i-1,j)+(1.0-b)*(circ_a(i,j+1)+circ_a(i,j-1)) &
           +b*(circ_aa(i,j+1)+circ_aa(i,j-1)))
       if(D <= 0.) goto 10
       DDD=0.5*(circ_a(i+1,j)+circ_a(i-1,j))

       f1=0.125*dtx*S*Um*(circ_a(i+2,j)-circ_a(i-2,j))
       f2=0.125*S*dtx*circ_a(i,j)
       f3=0.5*ocean_dt*ocf*Um
       f4=0.25*dtx*oc_gamma*(circ_p(i,j+1)-circ_p(i,j-1))/(rho/1000.0)
       f4=f4*alin
       f5=0.25*dtx*grav*(b*(circ_aa(i,j+1)-circ_aa(i,j-1))+(1.0-b)*&
            (circ_a(i,j+1)-circ_a(i,j-1)))

       ViXm=0.25*(circ_p(i+1,j-1)+circ_p(i+1,j+1)+circ_p(i-1,j+1) &
           +circ_p(i-1,j-1))
       VimOD=sqrt(circ_p(i,j)*circ_p(i,j)+ViXm*ViXm)
       VimOD=VimOD*51.44

       Ca=0.0012875
       if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD

       TSX=0.001*Ca*51.44*ViXm*VimOD
       TSY=0.001*Ca*51.44*circ_p(i,j)*VimOD
       TSX=TSX*alin
       TSY=TSY*alin

       f6=0.5*ocean_dt*TSY/(rho*D)

       uhm=Um
       vhm=b*circ_aa(i,j)+(1.0-b)*circ_a(i,j)

       f7=0.5*ocean_dt*oceanfric*sqrt(uhm*uhm+vhm*vhm)/D
       f8=0.125*ah*dtx*(circ_a(i+2,j)-2*circ_a(i,j)+circ_a(i-2,j)) /circdx
       f9=0.125*ah*dtx*(1.0-b)*(circ_a(i,j+2)-2*circ_a(i,j)+&
            circ_a(i,j-2))/circdx

       if(iSen1==1.and.j==ini) goto 50
       if(iSen2==1.and.j==ifi) goto 50

       goto 60

50     circ_aa(i,j)=(circ_a(i,j)*(1.0-(1.0-b)*f7)-f3-f4-f5+f6)/(1.0+b*f7)

       if(iSen1==1.and.j==ini)then
          ini=ini+2
          iSen1=0
          goto 30
       else
          ifi=ifi-2
          iSen2=0
          goto 30
       endif

60     g1=- b*(f2+0.125*dtx*ah/circdx)
       g2=1.0+b*(f7+0.25*dtx*ah/circdx)
       g3=b*(f2-0.125*dtx*ah/circdx)
       g4=circ_a(i,j)-f1-(1.0-b)*f2*(circ_a(i,j+2)-circ_a(i,j-2))-&
            f3-f4-f5+f6-f7*(1.0-b)*circ_a(i,j)+f8+f9

       if(j/=ini) goto 55
       U1(j-2)=0.0
       W1(j-2)=circ_aa(i,j-2)

55     X=g2-U1(j-2)*g1
       U1(j)=g3/X
       W1(j)=(g4-W1(j-2)*g1)/X
       if(j >= ifi) goto 15
       j=j+2
       goto 35

       ! Computes the VV values
15     do j=ifi,ini,-2
          circ_aa(i,j)=W1(j)-U1(j) *circ_aa(i,j+2)
       enddo

10     continue

    enddo

    return

  end subroutine calculateV
  ! =====================================================================================
  subroutine calculateU

    integer::j,i,kn,iSen1,iSen2,i1,i2,ideS1,ini,ifi

    real::S,b,dtx,Vm,D,DDD,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real::ViYm,VimOD,Ca,X,TSX,TSY,uhm,vhm

    real::U1(0:circm),W1(0:circm)

    U1=0.
    W1=0.

    dtx=ocean_dt/circdx
    S=0.
    b=oc_beta

    do kn=1,nn
       iSen1=0
       iSen2=0
       i1=in1(kn)+2
       i2=in2(kn)-2
       j=jn(kn)
       ideS1=in2(kn)-in1(kn)
       if(ideS1==2) goto 10
       ini=i1
       ifi=i2
       if(ip(i1)==i1) goto 20
       i1=i1-1
       ini=i1
       iSen1=1
20     if(ip(i2)==i2) goto 30
       i2=i2+1
       ifi=i2
       iSen2=1
30     i=ini

35     if(iSen1==0.and.iSen2==1.and.i==ini) i=ifi

       if(i<=2.or.i>circm-2) goto 10

       Vm=0.25*(b*(circ_aa(i+1,j+1)+circ_aa(i+1,j-1)+circ_aa(i-1,j-1)+&
            circ_aa(i-1,j+1))+(1-b)*(circ_a(i+1,j+1)+circ_a(i+1,j-1)+&
            circ_a(i-1,j-1)+circ_a(i-1,j+1)))

       D=0.5*(circ_a(i,j+1)+circ_a(i,j-1)+(1-b)*(circ_a(i+1,j)+circ_a(i-1,j)) &
           +b*(circ_aa(i+1,j)+circ_aa(i-1,j)))

       if(D <= 0.) goto 10
       DDD=0.5*(circ_a(i,j+1)+circ_a(i,j-1))

       f1=0.125*S*dtx*circ_a(i,j)
       f2=0.125*S*dtx*Vm*(circ_a(i,j+2)-circ_a(i,j-2))
       f3=0.5*ocean_dt*ocf*Vm

       ViYm=0.25*(circ_p(i+1,j-1)+circ_p(i+1,j+1)+circ_p(i-1,j+1) &
           +circ_p(i-1,j-1))
       VimOD=sqrt(circ_p(i,j)*circ_p(i,j)+ViYm*ViYm)
       VimOD=VimOD*51.44
       Ca=0.0012875

       if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD
       TSX=0.001*Ca*51.44*circ_p(i,j)*VimOD
       TSY=0.001*Ca*51.44*ViYm*VimOD
       TSX=TSX*alin
       TSY=TSY*alin

       f4=0.5*ocean_dt*TSX/(rho*D)

       uhm=b*circ_aa(i,j)+(1.0-b)*circ_a(i,j)
       vhm=Vm

       f5=0.5*ocean_dt*oceanfric*sqrt(uhm*uhm+vhm*vhm)/D
       f6=0.25*dtx*grav*(b*(circ_aa(i+1,j)-circ_aa(i-1,j))+(1.0-b)*&
            (circ_a(i+1,j)-circ_a(i-1,j)))
       f7=0.25*dtx*oc_gamma*(circ_p(i+1,j)-circ_p(i-1,j))/(rho/1000.0)
       f7=f7*alin
       f8=0.125*ah*dtx*(1.0-b)*(circ_a(i+2,j)-2.0*circ_a(i,j)+circ_a(i-2,j))/circdx
       f9=0.125*ah*dtx *(circ_a(i,j+2)-2.0*circ_a(i,j)+circ_a(i,j-2))/circdx

       if(iSen1==1.and.i==ini) goto 50
       if(iSen2==1.and.i==ifi) goto 50
       goto 60

50     circ_aa(i,j)=(circ_a(i,j)*(1.0-(1.0-b)*f5)+f3+f4-f6-f7)/(1.0+b*f5)
       if(iSen1==1.and.i==ini) then
          ini=ini+2
          iSen1=0
          goto 30
       else
          ifi=ifi-2
          iSen2=0
          goto 30
       endif

60     g1=-b*(f1+0.125*dtx*ah/circdx)
       g2=1.0+b*(f5+0.25*dtx*ah/circdx)
       g3=b*(f1-0.125*dtx*ah/circdx)
       g4=circ_a(i,j)-f2-(1.0-b)*f1*(circ_a(i+2,j)-circ_a(i-2,j))+&
            f3+f4-f6-f7-f5*(1.0-b)*circ_a(i,j)+f8+f9
       if(i/=ini) goto 55

       U1(i-2)=0.0
       W1(i-2)=circ_aa(i-2,j)

55     X=g2-U1(i-2)*g1
       U1(i)=g3/X
       W1(i)=(g4-W1(i-2)*g1)/X
       if(i >= ifi) goto 15

       ! Computes the UU values
       i=i+2
       goto 35
15     do i=ifi,ini,-2
          circ_aa(i,j)=W1(i)-U1(i)*circ_aa(i+2,j)
       enddo

10     continue
    enddo

    return

  end subroutine calculateU
  ! =====================================================================================
  ! =====================================================================================
end module ocean_circ
