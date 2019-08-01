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
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  OceanModel.f90
!
!    Description:  Ocean circulation model.
!
!        Version:  1.0
!        Created:  01/10/15 08:50:10
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module ocean_model

  use parallel
  use topology
  use parameters
  use external_forces
  use ocean_circ

  implicit none

contains

  ! =====================================================================================
  subroutine ocean_initialisation

    integer::step,i,j,ic,jc,p,m
    real(kind=8)::ochmax
    real(kind=8),dimension(nbnodes)::tempZ

    if(ocean_dx<dx) ocean_dx=real(dx)
    step=int(ocean_dx/dx)

    if(step>1)then
      nbox=int(nx/step)+2
      nboy=int(ny/step)+2
    else
      nbox=nx
      nboy=ny
    endif
    oc_beta=0.6666666
    oc_gamma=0.0
    circdx=ocean_dx*100.

    ! South west corner
    ocean_xo=real(minx)-ocean_dx
    ocean_yo=real(miny)-ocean_dx
    if(ocean_xo+(nbox-1)*ocean_dx>minx+nx*dx)nbox=nbox-1
    if(ocean_yo+(nboy-1)*ocean_dx>miny+ny*dx)nboy=nboy-1

    if(allocated(oceanX)) deallocate(oceanX)
    allocate(oceanX(nbox))
    do i=1,nbox
      oceanX(i)=(i-1)*ocean_dx+ocean_xo
    enddo

    if(allocated(oceanY)) deallocate(oceanY)
    allocate(oceanY(nboy))
    do i=1,nboy
      oceanY(i)=(i-1)*ocean_dx+ocean_yo
    enddo

    if(allocated(oceanZ)) deallocate(oceanZ)
    allocate(oceanZ(nbox,nboy))

    if(allocated(circU)) deallocate(circU)
    allocate(circU(nbox,nboy,season))
    if(allocated(circV)) deallocate(circV)
    allocate(circV(nbox,nboy,season))

    ! Read topography regular grid.
    p=0
    m=0
    ! nt*,nx,ny
    ! open(unit=17,file='newnodes.csv')
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
          tempZ(m)=rcoordZ(p)
          ! write(17,*)m,(i-1)*3000.,(j-1)*3000.,tempZ(m)
        endif
      enddo
    enddo
    ! close(17)
    ! stop
    ochmax=0.0
    j=1
    p=0
    do jc=2,nboy-1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbox-1
        p=p+1
        oceanZ(ic,jc)=real(tempZ(p))
        ochmax=max(gsea%actual_sea-oceanZ(ic,jc),ochmax)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo
    ! Update border
    oceanZ(2:nbox-1,1)=oceanZ(2:nbox-1,2)!+(oceanZ(2:nbox-1,2)-oceanZ(2:nbox-1,3))
    oceanZ(2:nbox-1,nboy)=oceanZ(2:nbox-1,nboy-1)!+(oceanZ(2:nbox-1,nboy-1)-oceanZ(2:nbox-1,nboy-2))
    oceanZ(1,2:nboy-1)=oceanZ(2,2:nboy-1)!+oceanZ(2,2:nboy-1)-oceanZ(3,2:nboy-1)
    oceanZ(nbox,2:nboy-1)=oceanZ(nbox-1,2:nboy-1)!+oceanZ(nbox-1,2:nboy-1)-oceanZ(nbox-2,2:nboy-1)

    ! Update corner
    oceanZ(1,1)=oceanZ(2,2)!+oceanZ(2,2)-oceanZ(3,3)
    oceanZ(1,nboy)=oceanZ(2,nboy-1)!+oceanZ(2,nboy-1)-oceanZ(3,nboy-2)
    oceanZ(nbox,1)=oceanZ(nbox-1,2)!+oceanZ(nbox-1,2)-oceanZ(nbox-2,3)
    oceanZ(nbox,nboy)=oceanZ(nbox-1,nboy-1)!+oceanZ(nbox-1,nboy-1)-oceanZ(nbox-2,nboy-2)

    ! Circulation computational grid
    circn=nbox+gap*2
    circm=nboy+gap*2
    if(ochmax==0.0)then
      ocean_dt=60.0
    else
      ocean_dt=real(Courant/(sqrt(9.81*ochmax)*sqrt(2/(ocean_dx**2))))
    endif
    ah=5.*circdx
    atf=100.*ocean_dt

    ocirc_steps=int(storm_dt*3600/ocean_dt)

    rho=1.0
    latitude=-latitude
    ome=2.*3.14/(24.*3600.)
    ocf=2.*ome*sin(latitude*3.14/180.)
    tiner=2.*3.14/abs(ocf)
    niner=tiner/ocean_dt

    ! Allocate forecasts
    forecast_param(1)=hindcast(1)%subgroup(1)%hs
    forecast_param(2)=hindcast(1)%subgroup(1)%per
    forecast_param(3)=hindcast(1)%subgroup(1)%dir
    forecast_param(4)=hindcast(1)%subgroup(1)%dd
    forecast_param(5)=hindcast(1)%subgroup(1)%wvel
    forecast_param(6)=hindcast(1)%subgroup(1)%wdir
    forecast_param(7)=forecast_param(5)*sin(forecast_param(6)*180./3.14)
    forecast_param(8)=forecast_param(5)*cos(forecast_param(6)* 180./3.14)

    if(pet_id==0)then
      if(allocated(circZ)) deallocate(circZ)
      if(allocated(circ_a)) deallocate(circ_a)
      if(allocated(circ_p)) deallocate(circ_p)
      if(allocated(circ_aa)) deallocate(circ_aa)
      allocate(circZ(circm,circn))
      allocate(circ_a(circm,circn),circ_aa(circm,circn),circ_p(circm,circn))
      if(allocated(p1)) deallocate(p1,p2,p3,p4,p5,p6,p7)
      allocate(p1(circn*circn,3),p2( circm*circm,3),p3(circn*circn,3),p4(circm*circm,3))
      allocate(p5(2*circm+2*circn,3),p6(circn*circn,3),p7(circn*circn,3))
    endif

    ! Arrays for interpolation
    if(allocated(x1_array)) deallocate(x1_array)
    if(allocated(y1_array)) deallocate(y1_array)
    allocate(x1_array(circn),y1_array(circm))

    if(allocated(xc_array)) deallocate(xc_array)
    if(allocated(yc_array)) deallocate(yc_array)
    allocate(xc_array(circn),yc_array(circm))
    ylen=0
    yclen=0
    y1_array=0.
    yc_array=0.
    do i=1,circm
      if(i==ip(i))then
        ylen=ylen+1
        y1_array(ylen)=(i-1)*real(dx)
      endif
      yclen=yclen+1
      yc_array(yclen)=(i-1)*real(dx)
    enddo

    xlen=0
    xclen=0
    x1_array=0.
    xc_array=0.
    do j=1,circn
      if(j==ip(j))then
        xlen=xlen+1
        x1_array(xlen)=(j-1)*real(dx)
      endif
      xclen=xclen+1
      xc_array(xclen)=(j-1)*real(dx)
    enddo

    if(allocated(o1)) deallocate(o1)
    if(allocated(o2)) deallocate(o2)
    allocate(o1(xlen,ylen),o2(xlen,ylen))

    if(.not.allocated(rcircU))then
      allocate(rcircU(nbox,nboy))
      allocate(rcircV(nbox,nboy))
    endif

    return

  end subroutine ocean_initialisation
  ! =====================================================================================
  subroutine ocean_reset_topography

    integer::step,i,j,ic,jc,p,m,ig,jg,li,lj,k
    real(kind=8),dimension(nbnodes)::tempZ

    if(ocean_dx<dx) ocean_dx=real(dx)
    step=int(ocean_dx/dx)

    if(pet_id/=0) return

    ! Read changes in topographic regular grid.
    p=0
    m=0
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
         tempZ(m)=rcoordZ(p)
        endif
      enddo
    enddo

    j=1
    do jc=2,nboy-1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbox-1
        oceanZ(ic,jc)=real(tempZ(p))
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    oceanZ(2:nbox-1,1)=oceanZ(2:nbox-1,2)!+(oceanZ(2:nbox-1,2)-oceanZ(2:nbox-1,3))
    oceanZ(2:nbox-1,nboy)=oceanZ(2:nbox-1,nboy-1)!+(oceanZ(2:nbox-1,nboy-1)-oceanZ(2:nbox-1,nboy-2))
    oceanZ(1,2:nboy-1)=oceanZ(2,2:nboy-1)!+oceanZ(2,2:nboy-1)-oceanZ(3,2:nboy-1)
    oceanZ(nbox,2:nboy-1)=oceanZ(nbox-1,2:nboy-1)!+oceanZ(nbox-1,2:nboy-1)-oceanZ(nbox-2,2:nboy-1)

    ! Update corner
    oceanZ(1,1)=oceanZ(2,2)!+oceanZ(2,2)-oceanZ(3,3)
    oceanZ(1,nboy)=oceanZ(2,nboy-1)!+oceanZ(2,nboy-1)-oceanZ(3,nboy-2)
    oceanZ(nbox,1)=oceanZ(nbox-1,2)!+oceanZ(nbox-1,2)-oceanZ(nbox-2,3)
    oceanZ(nbox,nboy)=oceanZ(nbox-1,nboy-1)!+oceanZ(nbox-1,nboy-1)-oceanZ(nbox-2,nboy-2)

    ! Update topography
    ig=0
    jg=0
    li=1
    lj=1
    do i=1,nbox
      do j=1,nboy
        oceanZ(i,j)=oceanZ(i,j)-real(gsea%actual_sea)
        if(ig==0.and.jg==0)then
          circZ(lj+gap,li+gap)=oceanZ(i,j)
          if(oceanZ(i,j)<=-circbase) circZ(lj+gap,li+gap)=-circbase
          lj=lj+1
        endif
        jg=jg+1
        if(jg==step) jg=0
      enddo
      ig=ig+1
      jg=0
      if(ig==step)then
        ig=0
        li=li+1
        lj=1
      endif
    enddo

    ! Update boundaries elevation
    do k=1,2
      do i=1,gap
        circZ(i,1:circn)=circZ(gap+1,1:circn)
        circZ(circm-i+1,1:circn)=circZ(circm-gap,1:circn)
      enddo
      do j=1,gap
        circZ(1:circm,j)=circZ(1:circm,gap+1)
        circZ(1:circm,circn-j+1)=circZ(1:circm,circn-gap)
      enddo
    enddo
    rcircU=0.
    rcircV=0.

    return

  end subroutine ocean_reset_topography
  ! =====================================================================================
  subroutine ocean_circulation_run(sgp)

    integer::i,j,sgp,p

    if(pet_id==0.and.circON==1)then
      ! Initialise computational arrays
      circ_a=0.0
      circ_p=0.0
      circ_aa=0.0

      ! Update topography
      p=1
      do i=1,circm
        do j=1,circn
          if(circZ(i,j)<0.)then
            circZ(i,j)=-circZ(i,j)*100.
          else
            circZ(i,j)=0.
          endif
          if(i==ip(i).and.j/=ip(j))then
            circ_a(i,j)=circZ(i,j)
          endif
        enddo
      enddo

      ! Specify the meteorological forces.
      do i=1,circm
        do j=1,circn
          ! Get wind velocity along x which needs to be converted in knots
          if(i==ip(i).and.j==ip(j))then
            circ_p(i,j)=forecast_param(7)*knots
          ! Get wind velocity along y which needs to be converted in knots
          elseif(i/=ip(i).and.j/=ip(j))then
            circ_p(i,j)=forecast_param(8)*knots
          endif
        enddo
      enddo

      ! Build boundaries
      call build_circulation_bounds

    endif
    call mpi_barrier(badlands_world,rc)

    call circulation_compute(sgp)

  end subroutine ocean_circulation_run
  ! =====================================================================================
  subroutine build_circulation_bounds

    logical::marine
    integer::i,j,st,minp1,minp2,halo

    halo=2

    ! Specify the positions i,j of the computational domain boundaries
    ! First compute boundaries for jn,in1,in2
    pt1=0
    p1=0
    minp1=circn
    do j=gap+1,circn-gap
       marine=.false.
       if(j==ip(j))then
          do i=gap+1,circm-gap
             if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and..not.marine.and.i/=circm-gap)then
                pt1=pt1+1
                p1(pt1,1)=j
                p1(pt1,2)=i
                marine=.true.
                minp1=min(minp1,j+1)
             elseif(marine.and.circZ(i,j)<=0.)then
                marine=.false.
                p1(pt1,3)=i-1
             elseif(marine.and.i==circm-gap)then
                p1(pt1,3)=i
             endif
          enddo
       endif
    enddo

    ! Then compute boundaries for im,jm1,jm2
    pt2=0
    p2=0
    minp2=circm
    do i=gap+1,circm-gap
       marine=.false.
       if(i/=ip(i))then
          do j=gap+1,circn-gap
             if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and..not.marine.and.j/=circn-gap)then
                pt2=pt2+1
                p2(pt2,1)=i
                p2(pt2,2)=j
                minp2=min(minp2,i+1)
                marine=.true.
             elseif(marine.and.circZ(i,j)<=0.)then
                marine=.false.
                p2(pt2,3)=j-1
             elseif(marine.and.j==circn-gap)then
                p2(pt2,3)=j
             endif
          enddo
       endif
    enddo

    ! Then compute boundaries for ina,inb,jnab
    pt3=0
    p3=0
    do j=gap+1,circn-gap
       marine=.false.
       do i=gap+1,circm-gap
          if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and..not.marine.and.i/=circm-gap)then
             pt3=pt3+1
             p3(pt3,1)=j
             p3(pt3,2)=i-1
             marine=.true.
          elseif(marine.and.circZ(i,j)<=0.)then
             marine=.false.
             p3(pt3,3)=i-1
          elseif(marine.and.i==circm-gap)then
             p3(pt3,3)=i-1
          endif
       enddo
    enddo

    ! Then compute boundaries for imab,jma,jmb
    pt4=0
    p4=0
    do i=gap+1,circm-gap
       marine=.false.
       do j=gap+1,circn-gap
          if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and..not.marine.and.j/=circn-gap)then
             pt4=pt4+1
             p4(pt4,1)=i
             p4(pt4,2)=j-1
             marine=.true.
          elseif(marine.and.circZ(i,j)<=0.)then
             marine=.false.
             p4(pt4,3)=j
          elseif(marine.and.j==circn-gap)then
             p4(pt4,3)=j
          endif
       enddo
    enddo

    ! Then compute boundaries for ie,je,ije
    ! Lower boundary
    pt5=0
    p5=0
    st=gap+1
    do j=gap+1,circn-gap
       if(circZ(st,j)>0.and.circZ(st,j)<circbase*100..and.j==ip(j))then
          pt5=pt5+1
          p5(pt5,1)=st
          p5(pt5,2)=j
          p5(pt5,3)=-1
       endif
    enddo
    ! Upper bondary
    st=circm-gap
    do j=gap+1,circn-gap
       if(circZ(st,j)>0.and.circZ(st,j)<circbase*100..and.j==ip(j))then
          pt5=pt5+1
          p5(pt5,1)=st
          p5(pt5,2)=j
          p5(pt5,3)=1
       endif
    enddo
    ! Left bondary
    st=gap+1
    do i=gap+1,circm-gap
       if(circZ(i,st)>0.and.circZ(i,st)<circbase*100..and.i==ip(i))then
          pt5=pt5+1
          p5(pt5,1)=i
          p5(pt5,2)=st
          p5(pt5,3)=-2
       endif
    enddo
    ! Right bondary
    st=circn-gap
    do i=gap+1,circm-gap
       if(circZ(i,st)>0.and.circZ(i,st)<circbase*100..and.i==ip(i))then
          pt5=pt5+1
          p5(pt5,1)=i
          p5(pt5,2)=st
          p5(pt5,3)=2
       endif
    enddo
    ! In case there is no open boundary
    if(pt5==0)then
       pt5=1
       p5(1,1)=1
       p5(1,2)=1
       p5(1,3)=0
    endif
    ! Then compute boundaries for inn1,inn2,jnn
    pt6=0
    p6=0
    if(minp1==0)print*,'Problem defining jnn'
    minp1=gap+1
    minp2=gap+1
    do j=minp1,circn-gap
       marine=.false.
       if(j/=ip(j))then
          do i=gap+1,circm-gap
             if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and. &
              .not.marine.and.i/=circm-gap)then
                pt6=pt6+1
                p6(pt6,1)=j
                p6(pt6,2)=i
                marine=.true.
             elseif(marine.and.circZ(i,j)<=0.)then
                marine=.false.
                p6(pt6,3)=i-1
             elseif(marine.and.i==circm-gap)then
                p6(pt6,3)=i
             endif
          enddo
       endif
    enddo
    ! Then compute boundaries for imm,jmm1,jmm2
    pt7=0
    p7=0
    if(minp2==0)print*,'Problem defining imm'
    do i=minp2,circm-gap
      marine=.false.
      if(i==ip(i))then
        do j=gap+1,circn-gap
          if(circZ(i,j)>0.and.circZ(i,j)<circbase*100..and..not.marine.and.j/=circn-gap)then
            pt7=pt7+1
            p7(pt7,1)=i
            p7(pt7,2)=j
            marine=.true.
          elseif(marine.and.circZ(i,j)<=0.)then
            marine=.false.
            p7(pt7,3)=j-1
          elseif(marine.and.j==circn-gap)then
            p7(pt7,3)=j
          endif
        enddo
      endif
    enddo

    ! Allocate computational grid boundaries
    nn=pt1
    if(allocated(jn)) deallocate(jn)
    if(allocated(in1)) deallocate(in1)
    if(allocated(in2)) deallocate(in2)
    allocate(in1(nn),in2(nn),jn(nn))
    in1=0
    in2=0
    jn=0
    jn(1:nn)=p1(1:nn,1)
    in1(1:nn)=p1(1:nn,2)
    in2(1:nn)=p1(1:nn,3)

    mm=pt2
    if(allocated(im)) deallocate(im)
    if(allocated(jm1)) deallocate(jm1)
    if(allocated(jm2)) deallocate(jm2)
    allocate(im(mm),jm1(mm),jm2(mm))
    im=0
    jm1=0
    jm2=0
    im(1:mm)=p2(1:mm,1)
    jm1(1:mm)=p2(1:mm,2)
    jm2(1:mm)=p2(1:mm,3)

    nnn=pt3
    if(allocated(ina)) deallocate(ina)
    if(allocated(inb)) deallocate(inb)
    if(allocated(jnab)) deallocate(jnab)
    allocate(ina(nnn),inb(nnn),jnab(nnn))
    ina=0
    inb=0
    jnab=0
    jnab(1:nnn)=p3(1:nnn,1)
    ina(1:nnn)=p3(1:nnn,2)
    inb(1:nnn)=p3(1:nnn,3)

    mmm=pt4
    if(allocated(jmb)) deallocate(jmb)
    if(allocated(jma)) deallocate(jma)
    if(allocated(imab)) deallocate(imab)
    allocate(imab(mmm),jma(mmm),jmb(mmm))
    imab=0
    jma=0
    jmb=0
    imab(1:mmm)=p4(1:mmm,1)
    jma(1:mmm)=p4(1:mmm,2)
    jmb(1:mmm)=p4(1:mmm,3)

    ne=pt5
    if(allocated(ie)) deallocate(ie)
    if(allocated(je)) deallocate(je)
    if(allocated(ije)) deallocate(ije)
    allocate(ie(ne),je(ne),ije(ne))
    ie=0
    je=0
    ije=0
    ie(1:ne)=p5(1:ne,1)
    je(1:ne)=p5(1:ne,2)
    ije(1:ne)=p5(1:ne,3)

    nnnn=pt6
    if(allocated(jnn)) deallocate(jnn)
    if(allocated(inn1)) deallocate(inn1)
    if(allocated(inn2)) deallocate(inn2)
    allocate(inn1(nnnn),inn2(nnnn),jnn(nnnn))
    inn1=0
    inn2=0
    jnn=0
    jnn(1:nnnn)=p6(1:nnnn,1)
    inn1(1:nnnn)=p6(1:nnnn,2)
    inn2(1:nnnn)=p6(1:nnnn,3)

    mmmm=pt7
    if(allocated(imm)) deallocate(imm)
    if(allocated(jmm1)) deallocate(jmm1)
    if(allocated(jmm2)) deallocate(jmm2)
    allocate(imm(mmmm),jmm1(mmmm), jmm2(mmmm))
    imm=0
    jmm1=0
    jmm2=0
    imm(1:mmmm)=p7(1:mmmm,1)
    jmm1(1:mmmm)=p7(1:mmmm,2)
    jmm2(1:mmmm)=p7(1:mmmm,3)

    return

  end subroutine build_circulation_bounds
  ! =====================================================================================
end module ocean_model
