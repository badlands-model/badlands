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
!       Filename:  IceModel.f90
!
!    Description:  Shallow Ice Approximation model.
!
!        Version:  1.0
!        Created:  01/07/15 10:14:30
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module ice_model

  use parallel
  use topology
  use parameters
  use gaussian_filter
  use external_forces

  implicit none

  integer::ice_stp

  real(kind=8)::gamma,gammas,ice_Hslide

  ! Glen's flow parameter
  integer,parameter::ice_glen=3

  ! Mass balance regularization parameter
  real(kind=8),parameter::ice_epsl=0.1
  ! Ice density
  real(kind=8),parameter::ice_dens=910.0

contains

  ! =====================================================================================
  subroutine ice_initialisation

    ! Initialise number of steps to reach steady step
    ice_stp=int(ice_Tstep/ice_dt)

    ! Define some computation coefficients
    gamma=2./(ice_glen+2)*ice_deform*(ice_dens*9.81)**ice_glen
    gammas=(ice_slide*ice_dens*9.81)**ice_glen

    ! Ice grid construction
    call  ICE_grid

    ! Smooth topography using Gaussian filter
    call gaussian_operator(iceZ,.true.,iceZb)
    iceH=0.0

    ! Update border
    iceZb(2:nbix-1,1)=iceZb(2:nbix-1,2)+(iceZb(2:nbix-1,2)-iceZb(2:nbix-1,3))
    iceZb(2:nbix-1,nbiy)=iceZb(2:nbix-1,nbiy-1)+(iceZb(2:nbix-1,nbiy-1)-iceZb(2:nbix-1,nbiy-2))
    iceZb(1,2:nbiy-1)=iceZb(2,2:nbiy-1)+iceZb(2,2:nbiy-1)-iceZb(3,2:nbiy-1)
    iceZb(nbix,2:nbiy-1)=iceZb(nbix-1,2:nbiy-1)+iceZb(nbix-1,2:nbiy-1)-iceZb(nbix-2,2:nbiy-1)

    ! Update corner
    iceZb(1,1)=iceZb(2,2)+iceZb(2,2)-iceZb(3,3)
    iceZb(1,nbiy)=iceZb(2,nbiy-1)+iceZb(2,nbiy-1)-iceZb(3,nbiy-2)
    iceZb(nbix,1)=iceZb(nbix-1,2)+iceZb(nbix-1,2)-iceZb(nbix-2,3)
    iceZb(nbix,nbiy)=iceZb(nbix-1,nbiy-1)+iceZb(nbix-1,nbiy-1)-iceZb(nbix-2,nbiy-2)

  end subroutine ice_initialisation
  ! =====================================================================================
  subroutine ice_reset_topography

    integer::step,i,j,ic,jc,p,m

    real(kind=8),dimension(nbnodes)::tempZ

    step=int(ice_dx/dx)

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
    do jc=2,nbiy-1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbix-1
        iceZ(ic,jc)=tempZ(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    iceZ(2:nbix-1,1)=iceZ(2:nbix-1,2)+(iceZ(2:nbix-1,2)-iceZ(2:nbix-1,3))
    iceZ(2:nbix-1,nbiy)=iceZ(2:nbix-1,nbiy-1)+(iceZ(2:nbix-1,nbiy-1)-iceZ(2:nbix-1,nbiy-2))
    iceZ(1,2:nbiy-1)=iceZ(2,2:nbiy-1)+iceZ(2,2:nbiy-1)-iceZ(3,2:nbiy-1)
    iceZ(nbix,2:nbiy-1)=iceZ(nbix-1,2:nbiy-1)+iceZ(nbix-1,2:nbiy-1)-iceZ(nbix-2,2:nbiy-1)

    ! Update corner
    iceZ(1,1)=iceZ(2,2)+iceZ(2,2)-iceZ(3,3)
    iceZ(1,nbiy)=iceZ(2,nbiy-1)+iceZ(2,nbiy-1)-iceZ(3,nbiy-2)
    iceZ(nbix,1)=iceZ(nbix-1,2)+iceZ(nbix-1,2)-iceZ(nbix-2,3)
    iceZ(nbix,nbiy)=iceZ(nbix-1,nbiy-1)+iceZ(nbix-1,nbiy-1)-iceZ(nbix-2,nbiy-2)

    ! Smooth topography using Gaussian filter
    call gaussian_operator(iceZ,.true.,iceZb)
    iceH=0.0

    ! Update border
    iceZb(2:nbix-1,1)=iceZb(2:nbix-1,2)+(iceZb(2:nbix-1,2)-iceZb(2:nbix-1,3))
    iceZb(2:nbix-1,nbiy)=iceZb(2:nbix-1,nbiy-1)+(iceZb(2:nbix-1,nbiy-1)-iceZb(2:nbix-1,nbiy-2))
    iceZb(1,2:nbiy-1)=iceZb(2,2:nbiy-1)+iceZb(2,2:nbiy-1)-iceZb(3,2:nbiy-1)
    iceZb(nbix,2:nbiy-1)=iceZb(nbix-1,2:nbiy-1)+iceZb(nbix-1,2:nbiy-1)-iceZb(nbix-2,2:nbiy-1)

    ! Update corner
    iceZb(1,1)=iceZb(2,2)+iceZb(2,2)-iceZb(3,3)
    iceZb(1,nbiy)=iceZb(2,nbiy-1)+iceZb(2,nbiy-1)-iceZb(3,nbiy-2)
    iceZb(nbix,1)=iceZb(nbix-1,2)+iceZb(nbix-1,2)-iceZb(nbix-2,3)
    iceZb(nbix,nbiy)=iceZb(nbix-1,nbiy-1)+iceZb(nbix-1,nbiy-1)-iceZb(nbix-2,nbiy-2)

  end subroutine ice_reset_topography
  ! =====================================================================================
  subroutine ice_SIA_run

    integer::i,j,p

    real(kind=8)::nldiff(4),a,b

    real(kind=8),dimension(nbix-2)::a1,b1,c1,d1,sol1
    real(kind=8),dimension(nbiy-2)::a2,b2,c2,d2,sol2

    real(kind=8),dimension(nbix,nbiy)::halfH
    real(kind=8),dimension(nbix*nbiy)::iceVU,iceVZ

    iceVU=0.
    iceVZ=0.

    if(pet_id==0)then

      ! Update Equilibrium Line Altitude for the considered time step
      if(gela%ela) call ELA_flux
      ! Update the sliding altitude accordingly
      ice_Hslide=ice_Zsld+gela%actual_ela

      ! Update top surface elevation
      iceZ=iceZb+iceH

  	  do p=1,ice_stp
        ! ADI: step one
        halfH=0.
        do j=2,nbiy-1
          do i=2,nbix-1
            nldiff=0.
            call compute_nonlinear_diffusion(i,j,nldiff)
            a1(i-1)=-ice_dt*(nldiff(1)+nldiff(2))/(4.0*ice_dx**2.)
            c1(i-1)=-ice_dt*(nldiff(3)+nldiff(4))/(4.0*ice_dx**2.)
            b1(i-1)=1.-(a1(i-1)+c1(i-1))
            a=-ice_dt*(nldiff(1)+nldiff(3))/(4.0*ice_dx**2.)
            b=-ice_dt*(nldiff(2)+nldiff(4))/(4.0*ice_dx**2.)
            d1(i-1)=iceZ(i,j)+surface_mass_balance(i,j)*ice_dt/2-( &
              a1(i-1)*iceZb(i-1,j)+b1(i-1)*iceZb(i,j)+c1(i-1)*iceZb(i+1,j))-( &
              a*iceZ(i,j-1)-(a+b)*iceZ(i,j)+b*iceZ(i,j+1))
            if(i==2) d1(i-1)=d1(i-1)-a1(i-1)*iceH(i-1,j)
            if(i==nbix-1) d1(i-1)=d1(i-1)-c1(i-1)*iceH(i+1,j)
          enddo
          call solve_tridiag(a1,b1,c1,d1,sol1,nbix-2)
          do i=2,nbix-1
            if(sol1(i-1)>0.)then
              halfH(i,j)=sol1(i-1)
            else
              halfH(i,j)=0.
            endif
          enddo
        enddo

        ! Update ice thickness boundaries
        call bc_zeroIce(halfH)
        ! Update ice thickness from first operator splitting time
        iceH=halfH
        ! Update top surface elevation consequently
        iceZ=iceZb+iceH

        ! ADI: step two
        halfH=0.
        do i=2,nbix-1
          do j=2,nbiy-1
            nldiff=0.
            call compute_nonlinear_diffusion(i,j,nldiff)
            a2(j-1)=-ice_dt*(nldiff(3)+nldiff(1))/(4.0*ice_dx**2.)
            c2(j-1)=-ice_dt*(nldiff(4)+nldiff(2))/(4.0*ice_dx**2.)
            b2(j-1)=1.-(a2(j-1)+c2(j-1))
            a=-ice_dt*(nldiff(1)+nldiff(2))/(4.0*ice_dx**2.)
            b=-ice_dt*(nldiff(3)+nldiff(4))/(4.0*ice_dx**2.)
            d2(j-1)=iceZ(i,j)+surface_mass_balance(i,j)*ice_dt/2-( &
              a2(j-1)*iceZb(i,j-1)+b2(j-1)*iceZb(i,j)+c2(j-1)*iceZb(i,j+1))-( &
              a*iceZ(i-1,j)-(a+b)*iceZ(i,j)+b*iceZ(i+1,j))
            if(j==2) d2(j-1)=d2(j-1)-a2(j-1)*iceH(i,j-1)
            if(j==nbiy-1) d2(j-1)=d2(j-1)-c2(j-1)*iceH(i,j+1)
          enddo
          call solve_tridiag(a2,b2,c2,d2,sol2,nbiy-2)
          do j=2,nbiy-1
            if(sol2(j-1)>0.)then
              halfH(i,j)=sol2(j-1)
            else
              halfH(i,j)=0.
            endif
          enddo
        enddo

        ! Update ice thickness boundaries
        call bc_zeroIce(halfH)
        ! Update ice thickness from second operator splitting time
        iceH=halfH
        ! Update top surface elevation consequently
        iceZ=iceZb+iceH

        iceU=0.
        do i=2,nbix-1
          do j=2,nbiy-1
            if(iceH(i,j)>0.0) iceU(i,j)=ice_sliding_velocity(i,j)
          enddo
        enddo

  	  enddo

      p=1
      do i=1,nbix
        do j=1,nbiy
          iceVU(p)=iceU(i,j)
          iceVZ(p)=iceZ(i,j)
          p=p+1
        enddo
      enddo

    endif

    call mpi_allreduce(mpi_in_place,iceVZ,nbix*nbiy,mpi_double_precision,mpi_max,badlands_world,rc)
    call mpi_allreduce(mpi_in_place,iceVU,nbix*nbiy,mpi_double_precision,mpi_max,badlands_world,rc)

    p=1
    do i=1,nbix
      do j=1,nbiy
        iceU(i,j)=iceVU(p)
        iceZ(i,j)=iceVZ(p)
        p=p+1
      enddo
    enddo


  end subroutine ice_SIA_run
  ! =====================================================================================
  subroutine compute_nonlinear_diffusion(i,j,nldiff)

    integer::i,j

    real(kind=8)::nnalpha12,ppalpha12,nnh12,pph12,nnkappa12,ppkappa12
    real(kind=8)::npalpha12,pnalpha12,nph12,pnh12,npkappa12,pnkappa12
    real(kind=8)::nngammas12,ppgammas12,npgammas12,pngammas12,nldiff(4)

    nnalpha12=((iceZ(i,j)-iceZ(i-1,j)+iceZ(i,j-1)-iceZ(i-1,j-1))/(2*ice_dx))**2
    nnalpha12=nnalpha12+((iceZ(i,j)-iceZ(i,j-1)+iceZ(i-1,j)-iceZ(i-1,j-1))/(2*ice_dx))**2
    nnalpha12=sqrt(nnalpha12)
    npalpha12=((iceZ(i,j)-iceZ(i-1,j)+iceZ(i,j+1)-iceZ(i-1,j+1))/(2*ice_dx))**2
    npalpha12=npalpha12+((iceZ(i,j)-iceZ(i,j+1)+iceZ(i-1,j)-iceZ(i-1,j+1))/(2*ice_dx))**2
    npalpha12=sqrt(npalpha12)
    pnalpha12=((iceZ(i,j)-iceZ(i+1,j)+iceZ(i,j-1)-iceZ(i+1,j-1))/(2*ice_dx))**2
    pnalpha12=pnalpha12+((iceZ(i,j)-iceZ(i,j-1)+iceZ(i+1,j)-iceZ(i+1,j-1))/(2*ice_dx))**2
    pnalpha12=sqrt(pnalpha12)
    ppalpha12=((iceZ(i,j)-iceZ(i+1,j)+iceZ(i,j+1)-iceZ(i+1,j+1))/(2*ice_dx))**2
    ppalpha12=ppalpha12+((iceZ(i,j)-iceZ(i,j+1)+iceZ(i+1,j)-iceZ(i+1,j+1))/(2*ice_dx))**2
    ppalpha12=sqrt(ppalpha12)

    nnh12=(iceH(i-1,j-1)+iceH(i,j-1)+iceH(i-1,j)+iceH(i,j))/4.
    nph12=(iceH(i-1,j+1)+iceH(i,j+1)+iceH(i-1,j)+iceH(i,j))/4.
    pnh12=(iceH(i+1,j-1)+iceH(i,j-1)+iceH(i+1,j)+iceH(i,j))/4.
    pph12=(iceH(i+1,j+1)+iceH(i,j+1)+iceH(i+1,j)+iceH(i,j))/4.

    nnkappa12=(nnh12**(ice_glen+1))*(nnalpha12**(ice_glen-1))
    npkappa12=(nph12**(ice_glen+1))*(npalpha12**(ice_glen-1))
    pnkappa12=(pnh12**(ice_glen+1))*(pnalpha12**(ice_glen-1))
    ppkappa12=(pph12**(ice_glen+1))*(ppalpha12**(ice_glen-1))

    nngammas12=0.0
    if(ice_Hslide>iceZb(i-1,j-1)) nngammas12=ice_Hslide-iceZb(i-1,j-1)
    if(ice_Hslide>iceZb(i,j-1)) nngammas12=nngammas12+ice_Hslide-iceZb(i,j-1)
    if(ice_Hslide>iceZb(i-1,j)) nngammas12=nngammas12+ice_Hslide-iceZb(i-1,j)
    if(ice_Hslide>iceZb(i,j)) nngammas12=nngammas12+ice_Hslide-iceZb(i,j)
    nngammas12=gammas*nngammas12/4.
    npgammas12=0.0
    if(ice_Hslide>iceZb(i-1,j+1)) npgammas12=ice_Hslide-iceZb(i-1,j+1)
    if(ice_Hslide>iceZb(i,j+1)) npgammas12=npgammas12+ice_Hslide-iceZb(i,j+1)
    if(ice_Hslide>iceZb(i-1,j)) npgammas12=npgammas12+ice_Hslide-iceZb(i-1,j)
    if(ice_Hslide>iceZb(i,j)) npgammas12=npgammas12+ice_Hslide-iceZb(i,j)
    npgammas12=gammas*npgammas12/4.
    pngammas12=0.0
    if(ice_Hslide>iceZb(i+1,j-1)) pngammas12=ice_Hslide-iceZb(i+1,j-1)
    if(ice_Hslide>iceZb(i,j-1)) pngammas12=pngammas12+ice_Hslide-iceZb(i,j-1)
    if(ice_Hslide>iceZb(i+1,j)) pngammas12=pngammas12+ice_Hslide-iceZb(i+1,j)
    if(ice_Hslide>iceZb(i,j)) pngammas12=pngammas12+ice_Hslide-iceZb(i,j)
    pngammas12=gammas*pngammas12/4.
    ppgammas12=0.0
    if(ice_Hslide>iceZb(i+1,j+1)) ppgammas12=ice_Hslide-iceZb(i+1,j+1)
    if(ice_Hslide>iceZb(i,j+1)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i,j+1)
    if(ice_Hslide>iceZb(i+1,j)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i+1,j)
    if(ice_Hslide>iceZb(i,j)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i,j)
    ppgammas12=gammas*ppgammas12/4.

    ! i-1/2 j-1/2
    nldiff(1)=(nnh12*gamma+nngammas12)*nnkappa12
    ! i-1/2 j+1/2
    nldiff(2)=(nph12*gamma+npgammas12)*npkappa12
    ! i+1/2 j-1/2
    nldiff(3)=(pnh12*gamma+pngammas12)*pnkappa12
    ! i+1/2 j+1/2
    nldiff(4)=(pph12*gamma+ppgammas12)*ppkappa12

  end subroutine compute_nonlinear_diffusion
  ! =====================================================================================
  function surface_mass_balance(i,j) result(mb)

    integer::i,j
    real(kind=8)::temp1,temp2
    real(kind=8)::mb

    temp1=0.5*(ice_m1+ice_m2)*(iceZ(i,j)-(gela%actual_ela+gsea%actual_sea))
    temp2=0.5*(ice_m1-ice_m2)*(iceZ(i,j)-(gela%actual_ela+gsea%actual_sea))

    mb=temp1-sqrt(temp2**2+ice_epsl)

  end function surface_mass_balance
  ! =====================================================================================
  subroutine solve_tridiag(a,b,c,d,x,n)

    ! a - sub-diagonal (means it is the diagonal below the main diagonal)
    ! b - the main diagonal
    ! c - sup-diagonal (means it is the diagonal above the main diagonal)
    ! d - right part
    ! x - the answer
    ! n - number of equations
    integer,intent(in)::n
    real(kind=8),dimension(n),intent(in)::a,b,c,d
    real(kind=8),dimension(n),intent(out)::x
    real(kind=8),dimension(n)::cp,dp
    real(kind=8)::m
    integer::i

    ! Initialize c-prime and d-prime
    cp(1)=c(1)/b(1)
    dp(1)=d(1)/b(1)

    ! Solve for vectors c-prime and d-prime
    do i=2,n
      m=b(i)-cp(i-1)*a(i)
      cp(i)=c(i)/m
      dp(i)=(d(i)-dp(i-1)*a(i))/m
    enddo

    ! Initialize x
    x(n)=dp(n)
    if(x(n)<0.)x(n)=0.

    ! Solve for x from the vectors c-prime and d-prime
    do i=n-1,1,-1
      x(i)=dp(i)-cp(i)*x(i+1)
      if(x(i)<0.)x(i)=0.
    enddo

  end subroutine solve_tridiag
  ! =====================================================================================
  function ice_sliding_velocity(i,j) result(us)

    integer::i,j

    real(kind=8)::nnalpha12,ppalpha12,nnh12,pph12,nnkappa12,ppkappa12
    real(kind=8)::npalpha12,pnalpha12,nph12,pnh12,npkappa12,pnkappa12
    real(kind=8)::nngammas12,ppgammas12,npgammas12,pngammas12

    real(kind=8)::us

    nnalpha12=((iceZ(i,j)-iceZ(i-1,j)+iceZ(i,j-1)-iceZ(i-1,j-1))/(2*ice_dx))**2
    nnalpha12=nnalpha12+((iceZ(i,j)-iceZ(i,j-1)+iceZ(i-1,j)-iceZ(i-1,j-1))/(2*ice_dx))**2
    nnalpha12=sqrt(nnalpha12)
    npalpha12=((iceZ(i,j)-iceZ(i-1,j)+iceZ(i,j+1)-iceZ(i-1,j+1))/(2*ice_dx))**2
    npalpha12=npalpha12+((iceZ(i,j)-iceZ(i,j+1)+iceZ(i-1,j)-iceZ(i-1,j+1))/(2*ice_dx))**2
    npalpha12=sqrt(npalpha12)
    pnalpha12=((iceZ(i,j)-iceZ(i+1,j)+iceZ(i,j-1)-iceZ(i+1,j-1))/(2*ice_dx))**2
    pnalpha12=pnalpha12+((iceZ(i,j)-iceZ(i,j-1)+iceZ(i+1,j)-iceZ(i+1,j-1))/(2*ice_dx))**2
    pnalpha12=sqrt(pnalpha12)
    ppalpha12=((iceZ(i,j)-iceZ(i+1,j)+iceZ(i,j+1)-iceZ(i+1,j+1))/(2*ice_dx))**2
    ppalpha12=ppalpha12+((iceZ(i,j)-iceZ(i,j+1)+iceZ(i+1,j)-iceZ(i+1,j+1))/(2*ice_dx))**2
    ppalpha12=sqrt(ppalpha12)

    nnh12=(iceH(i-1,j-1)+iceH(i,j-1)+iceH(i-1,j)+iceH(i,j))/4.
    nph12=(iceH(i-1,j+1)+iceH(i,j+1)+iceH(i-1,j)+iceH(i,j))/4.
    pnh12=(iceH(i+1,j-1)+iceH(i,j-1)+iceH(i+1,j)+iceH(i,j))/4.
    pph12=(iceH(i+1,j+1)+iceH(i,j+1)+iceH(i+1,j)+iceH(i,j))/4.

    nnkappa12=(nnh12**(ice_glen))*(nnalpha12**(ice_glen))
    npkappa12=(nph12**(ice_glen))*(npalpha12**(ice_glen))
    pnkappa12=(pnh12**(ice_glen))*(pnalpha12**(ice_glen))
    ppkappa12=(pph12**(ice_glen))*(ppalpha12**(ice_glen))

    nngammas12=0.0
    if(ice_Hslide>iceZb(i-1,j-1)) nngammas12=ice_Hslide-iceZb(i-1,j-1)
    if(ice_Hslide>iceZb(i,j-1)) nngammas12=nngammas12+ice_Hslide-iceZb(i,j-1)
    if(ice_Hslide>iceZb(i-1,j)) nngammas12=nngammas12+ice_Hslide-iceZb(i-1,j)
    if(ice_Hslide>iceZb(i,j)) nngammas12=nngammas12+ice_Hslide-iceZb(i,j)
    nngammas12=gammas*nngammas12/4.
    npgammas12=0.0
    if(ice_Hslide>iceZb(i-1,j+1)) npgammas12=ice_Hslide-iceZb(i-1,j+1)
    if(ice_Hslide>iceZb(i,j+1)) npgammas12=npgammas12+ice_Hslide-iceZb(i,j+1)
    if(ice_Hslide>iceZb(i-1,j)) npgammas12=npgammas12+ice_Hslide-iceZb(i-1,j)
    if(ice_Hslide>iceZb(i,j)) npgammas12=npgammas12+ice_Hslide-iceZb(i,j)
    npgammas12=gammas*npgammas12/4.
    pngammas12=0.0
    if(ice_Hslide>iceZb(i+1,j-1)) pngammas12=ice_Hslide-iceZb(i+1,j-1)
    if(ice_Hslide>iceZb(i,j-1)) pngammas12=pngammas12+ice_Hslide-iceZb(i,j-1)
    if(ice_Hslide>iceZb(i+1,j)) pngammas12=pngammas12+ice_Hslide-iceZb(i+1,j)
    if(ice_Hslide>iceZb(i,j)) pngammas12=pngammas12+ice_Hslide-iceZb(i,j)
    pngammas12=gammas*pngammas12/4.
    ppgammas12=0.0
    if(ice_Hslide>iceZb(i+1,j+1)) ppgammas12=ice_Hslide-iceZb(i+1,j+1)
    if(ice_Hslide>iceZb(i,j+1)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i,j+1)
    if(ice_Hslide>iceZb(i+1,j)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i+1,j)
    if(ice_Hslide>iceZb(i,j)) ppgammas12=ppgammas12+ice_Hslide-iceZb(i,j)
    ppgammas12=gammas*ppgammas12/4.

!     us=(nnh12*gamma*5./4.+nngammas12)*nnkappa12
!     us=us+(nph12*gamma*5./4.+npgammas12)*npkappa12
!     us=us+(pnh12*gamma*5./4.+pngammas12)*pnkappa12
!     us=us+(pph12*gamma*5./4.+ppgammas12)*ppkappa12
!     us=us/4.

    us=(nngammas12)*nnkappa12
    us=us+(npgammas12)*npkappa12
    us=us+(pngammas12)*pnkappa12
    us=us+(ppgammas12)*ppkappa12
    us=us/4.

  end function ice_sliding_velocity
  ! =====================================================================================

  subroutine bc_constantIce(bcArray)

    real(kind=8),intent(inout),dimension(nbix,nbiy)::bcArray

    ! Update border
    bcArray(2:nbix-1,1)=bcArray(2:nbix-1,2)
    bcArray(2:nbix-1,nbiy)=bcArray(2:nbix-1,nbiy-1)
    bcArray(1,2:nbiy-1)=bcArray(2,2:nbiy-1)
    bcArray(nbix,2:nbiy-1)=bcArray(nbix-1,2:nbiy-1)

    ! Update corner
    bcArray(1,1)=bcArray(2,2)
    bcArray(1,nbiy)=bcArray(2,nbiy-1)
    bcArray(nbix,1)=bcArray(nbix-1,2)
    bcArray(nbix,nbiy)=bcArray(nbix-1,nbiy-1)

  end subroutine bc_constantIce
  ! =====================================================================================
  subroutine bc_zeroIce(bcArray)

    real(kind=8),intent(inout),dimension(nbix,nbiy)::bcArray

    ! Update border
    bcArray(2:nbix-1,1)=0.
    bcArray(2:nbix-1,nbiy)=0.
    bcArray(1,2:nbiy-1)=0.
    bcArray(nbix,2:nbiy-1)=0.

    ! Update corner
    bcArray(1,1)=0.
    bcArray(1,nbiy)=0.
    bcArray(nbix,1)=0.
    bcArray(nbix,nbiy)=0.

  end subroutine bc_zeroIce
  ! =====================================================================================
  subroutine bc_wall

    ! Update border
    iceH(2:nbix-1,1)=0.
    iceH(2:nbix-1,nbiy)=0.
    iceH(1,2:nbiy-1)=0.
    iceH(nbix,2:nbiy-1)=0.

    ! Update corner
    iceH(1,1)=0.
    iceH(1,nbiy)=0.
    iceH(nbix,1)=0.
    iceH(nbix,nbiy)=0.

    ! Update border
    iceZb(2:nbix-1,1)=500.
    iceZb(2:nbix-1,nbiy)=500.
    iceZb(1,2:nbiy-1)=500.
    iceZb(nbix,2:nbiy-1)=500.

    ! Update corner
    iceZb(1,1)=500.
    iceZb(1,nbiy)=500.
    iceZb(nbix,1)=500.
    iceZb(nbix,nbiy)=500.

  end subroutine bc_wall
  ! =====================================================================================
end module ice_model
