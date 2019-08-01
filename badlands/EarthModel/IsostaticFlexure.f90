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
!       Filename:  IsostaticFlexure.f90
!
!    Description:  Flexural isostasy computation (Li et al. 2004 CompGeo) and compaction
!
!        Version:  1.0
!        Created:  11/07/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================

module isoflex

  use parallel
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none

  integer::numrow,numcol,n2row,n2col,n4row,n4col

  ! Fourth order schema parameters
  real(kind=8),parameter::q1=216.,q2=-18.,q3=-84.
  real(kind=8),parameter::q4=-3.,q5=792.
  real(kind=8),parameter::r1=144.,r2=18.,r3=-48.
  real(kind=8),parameter::r4=-6.,r5=1./240.

  real(kind=8),dimension(:,:),allocatable::w,wx,wy,ld,w1

contains

  ! =====================================================================================
  subroutine porosity_compaction

    integer::k,p

    real(kind=8)::phi,subs,pressure_lithos,mass,nh

    do k=1,dnodes
      pressure_lithos=0.
      if(spmZ(k)<gsea%actual_sea) pressure_lithos=9.81*sea_water_density*(gsea%actual_sea-spmZ(k))
      subs=0.
      do p=flex_lay-1,2,-1
        ! Mass due to water in pore space
        mass=mean_sediment_density*ulay_th(k,p+1)*(1-ulay_phi(k,p+1))
        mass=mass+sea_water_density*ulay_th(k,p+1)*ulay_phi(k,p+1)
        ! There is some sediment in this layer
        if(mass>0.)then
          pressure_lithos=pressure_lithos+mass*9.81
          ! Calculate new porosity with this pressure
          call get_porosity(pressure_lithos,ulay_phi(k,p),phi)
          if(phi>ulay_phi(k,p)) print*,'Issue computing compactional porosity'
          phi=min(phi,ulay_phi(k,p))
          if(phi<1.)then
            nh=ulay_th(k,p)*(1.0-ulay_phi(k,p))/(1.0-phi)
            nh=min(ulay_th(k,p),nh)
            ! Subsidence due to porosity change
            subs=subs+nh-ulay_th(k,p)
            ! Update layer thickness
            ulay_th(k,p)=nh
            ! Update porosity
            ulay_phi(k,p)=phi
          else
            ulay_th(k,p)=0.
            ulay_phi(k,p)=0.
          endif
        endif
      enddo
      ! Correct the topographic elevation due to compactional subsidence
      if(subs>0.0001) print*,'Problem when updating compactional subsidence'
      if(subs>0.) subs=0.
      spmZ(k)=spmZ(k)+subs
      sedthick(k)=sedthick(k)+subs
    enddo

  end subroutine porosity_compaction
  ! =====================================================================================
  subroutine get_porosity(Plith,in_phi,out_phi)

    integer:: p,fd
    real(kind=8)::in_phi,out_phi,Plith

    ! Go through the pressure field and find the fitting interval
    if(Plith<pressTable(1))then
      out_phi=poroTable(1)
    endif
    fd=0
    loop: do p=2,pressureFields
      if(Plith>=pressTable(p-1).and.Plith<pressTable(p))then
        fd=p
        out_phi=(Plith-pressTable(p-1))/(pressTable(p)-pressTable(p-1))
        out_phi=out_phi*(poroTable(p)-poroTable(p-1))+poroTable(p-1)
        exit loop
      endif
    enddo loop
    if(fd==0) out_phi=poroTable(pressureFields)

    out_phi=min(out_phi,in_phi)

  end subroutine get_porosity
  ! =====================================================================================
  subroutine update_Flex_array

    integer::step,i,j,ic,jc,p,m
    real(kind=8),dimension(nbnodes)::tempZ,tempSh

    step=int(flex_dx/dx)

    ! Read changes in topographic regular grid.
    p=0
    m=0
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
          tempZ(m)=rcoordZ(p)
          tempSh(m)=rsedload(p)
        endif
      enddo
    enddo

    j=1
    do jc=2,nbfy+1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbfx+1
        flexZ(ic,jc)=tempZ(p)
        flexSed(ic,jc)=tempSh(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    flexZ(2:nbfx+1,1)=flexZ(2:nbfx+1,2)
    flexZ(2:nbfx+1,nbfy+2)=flexZ(2:nbfx+1,nbfy+1)
    flexZ(1,2:nbfy+1)=flexZ(2,2:nbfy+1)
    flexZ(nbfx+2,2:nbfy+1)=flexZ(nbfx+1,2:nbfy+1)

    flexSed(2:nbfx+1,1)=flexSed(2:nbfx+1,2)
    flexSed(2:nbfx+1,nbfy+2)=flexSed(2:nbfx+1,nbfy+1)
    flexSed(1,2:nbfy+1)=flexSed(2,2:nbfy+1)
    flexSed(nbfx+2,2:nbfy+1)=flexSed(nbfx+1,2:nbfy+1)

    ! Update corner
    flexZ(1,1)=flexZ(2,2)
    flexZ(1,nbfy+2)=flexZ(2,nbfy+1)
    flexZ(nbfx+2,1)=flexZ(nbfx+1,2)
    flexZ(nbfx+2,nbfy+2)=flexZ(nbfx+1,nbfy+1)

    flexSed(1,1)=flexSed(2,2)
    flexSed(1,nbfy+2)=flexSed(2,nbfy+1)
    flexSed(nbfx+2,1)=flexSed(nbfx+1,2)
    flexSed(nbfx+2,nbfy+2)=flexSed(nbfx+1,nbfy+1)

  end subroutine update_Flex_array
  ! =====================================================================================
  subroutine built_initial_load

    integer::step,i,j,ic,jc,p,m

    real(kind=8)::tmp
    real(kind=8),dimension(nbnodes)::tempZ

    step=int(flex_dx/dx)

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
    do jc=2,nbfy+1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbfx+1
        flexZ(ic,jc)=tempZ(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border 
    flexZ(2:nbfx+1,1)=flexZ(2:nbfx+1,2)
    flexZ(2:nbfx+1,nbfy+2)=flexZ(2:nbfx+1,nbfy+1)
    flexZ(1,2:nbfy+1)=flexZ(2,2:nbfy+1)
    flexZ(nbfx+2,2:nbfy+1)=flexZ(nbfx+1,2:nbfy+1)

    ! Update corner
    flexZ(1,1)=flexZ(2,2)
    flexZ(1,nbfy+2)=flexZ(2,nbfy+1)
    flexZ(nbfx+2,1)=flexZ(nbfx+1,2)
    flexZ(nbfx+2,nbfy+2)=flexZ(nbfx+1,nbfy+1)

    if(allocated(sedloader))then
      p=0
      do j=1,nbfy+2
        do i=1,nbfx+2
          p=p+1
          prevload(i,j)=sedloader(p)
        enddo
      enddo
    else
      do j=1,nbfy+2
        do i=1,nbfx+2
          flexSed(i,j)=mean_sediment_density*100000.0*(1-poroTable(pressureFields))+100000.0*poroTable(pressureFields)*sea_water_density
          tmp=flexSed(i,j) !mean_sediment_density*flexSed(i,j) !*(1-flexPor(i,j))+flexPor(i,j)*flexSed(i,j)*sea_water_density
          if(flexZ(i,j)<gsea%actual_sea) tmp=tmp+(gsea%actual_sea-flexZ(i,j))*sea_water_density
          prevload(i,j)=tmp*cst1
        enddo
      enddo
    endif



  end subroutine built_initial_load
  ! =====================================================================================
  subroutine isostatic_flexure

    integer::i,j,n,m

    real(kind=8)::tmp,oldload,dtot,wtot,wdiff

    load=0.
    dtot=0.0
    do j=1,nbfy+2
      do i=1,nbfx+2
        if(j>1.and.j<nbfy+2.and.i>1.and.i<nbfx+2)then
          tmp=flexSed(i,j)
          if(flexZ(i,j)<gsea%actual_sea) tmp=tmp+(gsea%actual_sea-flexZ(i,j))*sea_water_density
          oldload=prevload(i,j)
          prevload(i,j)=cst1*tmp
          load(i,j)=prevload(i,j)-oldload
          dtot=dtot+load(i,j)
        endif
      enddo
      if(j>1.and.j<nbfy+2)then
        prevload(1,j)=prevload(2,j)
        prevload(nbfx+2,j)=prevload(nbfx+1,j)
        load(1,j)=load(2,j)
        load(nbfx+2,j)=load(nbfx+1,j)
      endif
    enddo

    ! Tiny added load, ignore isostatic changes
    if(abs(dtot)<cst1*mean_sediment_density)return

    prevload(1:nbfx+2,1)=prevload(1:nbfx+2,2)
    prevload(1:nbfx+2,nbfy+2)=prevload(1:nbfx+2,nbfy+1)
    load(1:nbfx+2,1)=load(1:nbfx+2,2)
    load(1:nbfx+2,nbfy+2)=load(1:nbfx+2,nbfy+1)

    ! Initialise flexure arrays
    numrow=nbfx+2
    numcol=nbfy+2
    if(.not.allocated(w)) allocate(w(numrow,numcol))
    if(.not.allocated(w1)) allocate(w1(numrow,numcol))
    if(.not.allocated(wx)) allocate(wx(numrow,numcol))
    if(.not.allocated(wy)) allocate(wy(numrow,numcol))
    if(.not.allocated(ld)) allocate(ld(numrow,numcol))

    ! Update flexure parameters
    w=0.0
    wx=0.0
    wy=0.0
    w1=0.0
    ld=load

    n2col=(numcol-2)/2
    n2row=(numrow-2)/2
    n4col=(numcol-2)/4
    n4row=(numrow-2)/4
    n2col=n2col*2+2
    n2row=n2row*2+2
    n4col=n4col*4+2
    n4row=n4row*4+2

    n=0
    ! Iterate until a solution is reached
    do
      n=n+1
      ! Fine grid first
      m=1
      if(flexorder==2)then
        call solve_flexure2(m,numrow-1,numcol-1,4)
      else
        call solve_flexure4(m,numrow-1,numcol-1,4)
      endif
      wtot=0.0
      wdiff=0.0
      do j=2,numcol-1
        do i=2,numrow-1
          wdiff=wdiff+abs(w1(i,j)-w(i,j))
          wtot=wtot+abs(w(i,j))
        enddo
      enddo
      w1=w
      if(wdiff<torb*wtot.or.n>100000)exit
      if(wtot==0.0)exit

      ! Coarser grid 2dx, full weighting operator
      m=2
      do j=2+m,n2col-m,m
        do i=2+m,n2row-m,m
          w(i,j)=0.0625*(w(i-1,j-1)+w(i+1,j-1)+w(i-1,j+1)+w(i+1,j+1)+ &
            2.0*(w(i,j-1)+w(i,j+1)+w(i-1,j)+w(i+1,j))+4.0*w(i,j))
          wx(i,j)=0.0625*(wx(i-1,j-1)+wx(i+1,j-1)+wx(i-1,j+1)+wx(i+1,j+1)+ &
            2.0*(wx(i,j-1)+wx(i,j+1)+wx(i-1,j)+wx(i+1,j))+4.0*wx(i,j))
          wy(i,j)=0.0625*(wy(i-1,j-1)+wy(i+1,j-1)+wy(i-1,j+1)+wy(i+1,j+1)+ &
            2.0*(wy(i,j-1)+wy(i,j+1)+wy(i-1,j)+wy(i+1,j))+4.0*wy(i,j))
        enddo
      enddo
      if(flexorder==2)then
        call solve_flexure2(m,n2row,n2col,4)
      else
        call solve_flexure4(m,n2row,n2col,4)
      endif
      ! Coarser grid 4dx, full weighting operator
      m=4
      do j=2+m,n4col-m,m
        do i=2+m,n4row-m,m
          w(i,j)=0.0625*(w(i-2,j-2)+w(i+2,j-2)+w(i-2,j+2)+w(i+2,j+2)+ &
            2.0*(w(i,j-2)+w(i,j+2)+w(i-2,j)+w(i+2,j))+4.0*w(i,j))
          wx(i,j)=0.0625*(wx(i-2,j-2)+wx(i+2,j-2)+wx(i-2,j+2)+wx(i+2,j+2)+ &
            2.0*(wx(i,j-2)+wx(i,j+2)+wx(i-2,j)+wx(i+2,j))+4.0*wx(i,j))
          wy(i,j)=0.0625*(wy(i-2,j-2)+wy(i+2,j-2)+wy(i-2,j+2)+wy(i+2,j+2)+ &
            2.0*(wy(i,j-2)+wy(i,j+2)+wy(i-2,j)+w(i+2,j))+4.0*wy(i,j))
        enddo
      enddo
      call boundary_flexure(m,n4row-m,n4col-m,w)
      call boundary_flexure(m,n4row-m,n4col-m,wx)
      call boundary_flexure(m,n4row-m,n4col-m,wy)
      if(flexorder==2)then
        call solve_flexure2(m,n4row,n4col,16)
      else
        call solve_flexure4(m,n4row,n4col,16)
      endif

      ! Interpolate to finer grid (2dx)
      do j=2,n4col-m,4
        do i=2,n4row-m,4
          w(i,j+2)=0.5*(w(i,j)+w(i,j+m))
          w(i+2,j)=0.5*(w(i,j)+w(i+m,j))
          w(i+2,j+2)=0.25*(w(i,j)+w(i,j+m)+w(i+m,j)+w(i+m,j+m))
          wx(i,j+2)=0.5*(wx(i,j)+wx(i,j+m))
          wx(i+2,j)=0.5*(wx(i,j)+wx(i+m,j))
          wx(i+2,j+2)=0.25*(wx(i,j)+wx(i,j+m)+wx(i+m,j)+wx(i+m,j+m))
          wy(i,j+2)=0.5*(wy(i,j)+wy(i,j+m))
          wy(i+2,j)=0.5*(wy(i,j)+wy(i+m,j))
          wy(i+2,j+2)=0.25*(wy(i,j)+wy(i,j+m)+wy(i+m,j)+wy(i+m,j+m))
        enddo
      enddo

      ! At boundaries of finer grid (2dx)
      if(n2col>n4col)then
        do i=2,n4row-m,4
          w(i,n2col)=w(i,n4col)
          w(i+2,n2col)=0.5*(w(i,n4col)+w(i+m,n4col))
          wx(i,n2col)=wx(i,n4col)
          wx(i+2,n2col)=0.5*(wx(i,n4col)+wx(i+m,n4col))
          wy(i,n2col)=wy(i,n4col)
          wy(i+2,n2col)=0.5*(wy(i,n4col)+wy(i+m,n4col))
        enddo
      endif

      if(n2row>n4row)then
        do j=2,n4col-m,4
          w(n2row,j)=w(n4row,j)
          w(n2row,j+2)=0.5*(w(n4row,j)+w(n4row,j+m))
          wx(n2row,j)=wx(n4row,j)
          wx(n2row,j+2)=0.5*(wx(n4row,j)+wx(n4row,j+m))
          wy(n2row,j)=wy(n4row,j)
          wy(n2row,j+2)=0.5*(wy(n4row,j)+wy(n4row,j+m))
        enddo
      endif
      w(n2row,n2col)=w(n4row,n4col)
      wx(n2row,n2col)=wx(n4row,n4col)
      wy(n2row,n2col)=wy(n4row,n4col)
      m=2
      if(flexorder==2)then
        call solve_flexure2(m,n2row,n2col,4)
      else
        call solve_flexure4(m,n2row,n2col,4)
      endif

      ! Finer grid dx
      do j=2,n2col-m,2
        do i=2,n2row-m,2
          w(i,j+1)=0.5*(w(i,j)+w(i,j+m))
          w(i+1,j)=0.5*(w(i,j)+w(i+m,j))
          w(i+1,j+1)=0.25*(w(i,j)+w(i,j+m)+w(i+m,j)+w(i+m,j+m))
          wx(i,j+1)=0.5*(wx(i,j)+wx(i,j+m))
          wx(i+1,j)=0.5*(wx(i,j)+wx(i+m,j))
          wx(i+1,j+1)=0.25*(wx(i,j)+wx(i,j+m)+wx(i+m,j)+wx(i+m,j+m))
          wy(i,j+1)=0.5*(wy(i,j)+wy(i,j+m))
          wy(i+1,j)=0.5*(wy(i,j)+wy(i+m,j))
          wy(i+1,j+1)=0.25*(wy(i,j)+wy(i,j+m)+wy(i+m,j)+wy(i+m,j+m))
        enddo
      enddo

      if(numcol-1>n2col)then
        do i=2,n2row-m,2
          w(i,numcol-1)=w(i,n2col)
          w(i+1,numcol-1)=0.5*(w(i,n2col)+w(i+m,n2col))
          wx(i,numcol-1)=wx(i,n2col)
          wx(i+1,numcol-1)=0.5*(wx(i,n2col)+wx(i+m,n2col))
          wy(i,numcol-1)=wy(i,n2col)
          wy(i+1,numcol-1)=0.5*(wy(i,n2col)+wy(i+m,n2col))
        enddo
      endif
      if(numrow-1>n2row)then
        do j=2,n2col-m,2
          w(numrow-1,j)=w(n2row,j)
          w(numrow-1,j+1)=0.5*(w(n2row,j)+w(n2row,j+m))
          wx(numrow-1,j)=wx(n2row,j)
          wx(numrow-1,j+1)=0.5*(wx(n2row,j)+wx(n2row,j+m))
          wy(numrow-1,j)=wy(n2row,j)
          wy(numrow-1,j+1)=0.5*(wy(n2row,j)+wy(n2row,j+m))
        enddo
      endif
      w(numrow-1,numcol-1)=w(n2row,n2col)
      wx(numrow-1,numcol-1)=wx(n2row,n2col)
      wy(numrow-1,numcol-1)=wy(n2row,n2col)
      do i=1,numrow
        w(i,1)=w(i,2)
        w(i,numcol)=w(i,numcol-1)
        wx(i,1)=wx(i,2)
        wx(i,numcol)=wx(i,numcol-1)
        wy(i,1)=wy(i,2)
        wy(i,numcol)=wy(i,numcol-1)
      enddo
      do j=1,numcol
        w(1,j)=w(2,j)
        w(numrow,j)=w(numrow-1,j)
        wx(1,j)=wx(2,j)
        wx(numrow,j)=wx(numrow-1,j)
        wy(1,j)=wy(2,j)
        wy(numrow,j)=wy(numrow-1,j)
      enddo
    enddo

    if(n>100000.and.pet_id==0) print*,'Isostacy did not converge to a solution.'
    flexDisp=w

  end subroutine isostatic_flexure
  ! =====================================================================================
  subroutine boundary_flexure(ks,nrw,ncl,temp)

    integer::ncl,nrw,i,j,ks,k1
    real(kind=8),dimension(numrow,numcol)::temp

    k1=2+ks
    do i=1,numrow
      do j=1,k1-1
        temp(i,j)=temp(i,k1)
      enddo
      do j=ncl+1,numcol
        temp(i,j)=temp(i,ncl)
      enddo
    enddo

    do j=1,numcol
      do i=1,k1-1
        temp(i,j)=temp(k1,j)
      enddo
      do i=nrw+1,numrow
        temp(i,j)=temp(nrw,j)
      enddo
    enddo

  end subroutine boundary_flexure
  ! =====================================================================================
  subroutine solve_flexure4(m,nrw,ncl,nloop)

    integer::m,ncl,nrw,nloop,i,j,ks,isw,jsw,ipass,n

    real(kind=8)::w3,wtot,wdiff,piv,dxm4,el1

    ks=0
    if(m>1)ks=m
    dxm4=(flex_dx**4)*m**4

    do j=2,ncl,m
      do i=2,nrw,m
        el1=flexZ(i,j)
        if(el1<=gsea%actual_sea)then
          ld(i,j)=(load(i,j)-cst2*w(i,j))*dxm4
        else
          ld(i,j)=(load(i,j)-cst3*w(i,j))*dxm4
        endif
      enddo
    enddo
    call boundary_flexure(0,nrw,ncl,ld)
    ! call boundary_flexure(ks,nrw-ks,ncl-ks,ld)

    n=0
    do
      n=n+1
      wtot=0.0
      wdiff=0.0
      jsw=1
      do ipass=1,2
        isw=jsw
        do j=2+ks,ncl-ks,m
          do i=isw+1+ks,nrw-ks,2*m
            el1=flexZ(i,j)
            if(el1<=gsea%actual_sea)then
              piv=q5+11.0*cst2*dxm4
            else
              piv=q5+11.0*cst3*dxm4
            endif

            wx(i,j)=(r1*(w(i+m,j)-w(i-m,j)) &
              +r2*(w(i+m,j+m)-w(i-m,j+m)+w(i+m,j-m)-w(i-m,j-m))  &
              +r3*(wx(i+m,j)+wx(i-m,j)) &
              +r4*(wx(i+m,j+m)+wx(i-m,j+m)+wx(i+m,j-m)+wx(i-m,j-m)) &
              +r4*(wy(i+m,j+m)-wy(i-m,j+m)-wy(i+m,j-m)+wy(i-m,j-m)) &
              +(ld(i+m,j)-ld(i-m,j)))*r5

            wy(i,j)=(r1*(w(i,j+m)-w(i,j-m)) &
              +r2*(w(i+m,j+m)+w(i-m,j+m)-w(i+m,j-m)-w(i-m,j-m))  &
              +r3*(wy(i,j+m)+wy(i,j-m)) &
              +r4*(wy(i+m,j+m)+wy(i-m,j+m)+wy(i+m,j-m)+wy(i-m,j-m)) &
              +r4*(wx(i+m,j+m)-wx(i-m,j+m)-wx(i+m,j-m)+wx(i-m,j-m)) &
              +(ld(i,j+m)-ld(i,j-m)))*r5

            w3=(q1*(w(i+m,j)+w(i-m,j)+w(i,j+m)+w(i,j-m)) &
              +q2*(w(i+m,j+m)+w(i-m,j+m)+w(i+m,j-m)+w(i-m,j-m)) &
              +q3*(wx(i+m,j)-wx(i-m,j)+wy(i,j+m)-wy(i,j-m)) &
              +q4*(wx(i+m,j+m)-wx(i-m,j+m)+wx(i+m,j-m)-wx(i-m,j-m)) &
              +q4*(wy(i+m,j+m)+wy(i-m,j+m)-wy(i+m,j-m)-wy(i-m,j-m)) &
              +(ld(i+m,j)+ld(i-m,j)+ld(i,j+m)+ld(i,j-m)) &
              +11.0*load(i,j)*dxm4)/piv

            if(el1<=gsea%actual_sea)then
              ld(i,j)=(load(i,j)-cst2*w3)*dxm4
            else
              ld(i,j)=(load(i,j)-cst3*w3)*dxm4
            endif

            wdiff=wdiff+abs(w3-w(i,j))
            wtot=wtot+abs(w3)
            w(i,j)=w3
          enddo
          isw=m+2-isw
        enddo
        jsw=m+2-jsw
        call boundary_flexure(ks,nrw-ks,ncl-ks,w)
        call boundary_flexure(ks,nrw-ks,ncl-ks,wx)
        call boundary_flexure(ks,nrw-ks,ncl-ks,wy)
        call boundary_flexure(ks,nrw-ks,ncl-ks,ld)
      enddo

      if(wdiff<1.0E-10*wtot.or.n>nloop)exit
      if(wtot==0.0)exit
    enddo

  end subroutine solve_flexure4
  ! =====================================================================================
  subroutine solve_flexure2(m,nrw,ncl,nloop)

    integer::m,ncl,nrw,nloop,i,j,ks,im,jm,im2,jm2,n,ip,jp,ip2,jp2

    real(kind=8)::w3,wtot,wdiff,piv,dxm4,el1

    ks=0
    if(m>1)ks=m
    dxm4=(flex_dx**4)*m**4

    do j=2,ncl,m
      do i=2,nrw,m
        el1=flexZ(i,j)
        if(el1<=gsea%actual_sea)then
          ld(i,j)=(load(i,j)-cst2*w(i,j))*dxm4
        else
          ld(i,j)=(load(i,j)-cst3*w(i,j))*dxm4
        endif
      enddo
    enddo
    call boundary_flexure(ks,nrw,ncl,ld)

    n=0
    do
      n=n+1
      wtot=0.0
      wdiff=0.0
      do j=2+ks,ncl-ks,m
        do i=2+ks,nrw-ks,m
          el1=flexZ(i,j)
          if(el1<=gsea%actual_sea)then
            piv=20.0+11.0/15.0*cst2*dxm4
          else
            piv=20.0+11.0/15.0*cst3*dxm4
          endif

          im=i-m
        	im2=i-2*m
        	jm=j-m
        	jm2=j-2*m
        	if(i-m<1)im=1
        	if(j-m<1)jm=1
        	if(i-2*m<1)im2=1
        	if(j-2*m<1)jm2=1
        	ip=i+m
        	jp=j+m
        	ip2=i+2*m
        	jp2=j+2*m
        	if(i+m>nrw)ip=nrw
        	if(j+m>ncl)jp=ncl
        	if(i+2*m>nrw)ip2=nrw
        	if(j+2*m>ncl)jp2=ncl

          w3=(8.0_8*(w(ip,j)+w(im,j)+w(i,jp)+w(i,jm))	&
            -2.0_8*(w(ip,jp)+w(im,jp)+w(ip,jm)+w(im,jm)) &
            -(w(ip2,j)+w(im2,j)+w(i,jp2)+w(i,jm2))+ &
          	(ld(ip,j)+ld(im,j)+ld(i,jp)+ld(i,jm))/15.0_8 &
            +11.0_8*load(i,j)*dxm4/15.0_8)/piv

          if(el1<=gsea%actual_sea)then
            ld(i,j)=(load(i,j)-cst2*w3)*dxm4
          else
            ld(i,j)=(load(i,j)-cst3*w3)*dxm4
          endif
          wdiff=wdiff+abs(w3-w(i,j))
          wtot=wtot+abs(w3)
          w(i,j)=w3
        enddo
      enddo

      if(wdiff<1.0E-10*wtot.or.n>nloop)exit
      if(wtot==0.0)exit
      call boundary_flexure(ks,nrw,ncl,w)
      call boundary_flexure(ks,nrw,ncl,ld)
    enddo

  end subroutine solve_flexure2
  ! =====================================================================================
end module isoflex
! =====================================================================================
