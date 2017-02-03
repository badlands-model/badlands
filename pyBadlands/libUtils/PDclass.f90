!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm class
module pdclass

    implicit none

    integer :: dnodes

    ! Define the data-structure to hold the data
    integer :: s1 = 0
    integer :: s2 = 0
    integer, allocatable, dimension(:) :: data1
    integer, allocatable, dimension(:) :: data2

    ! Set the size of allocated memory blocks
    integer :: block_size

    ! Set neighbourhood array
    integer,allocatable, dimension(:) :: maxngb
    integer,allocatable, dimension(:,:) :: neighbours

    ! Set edge distance array
    real(kind=8),allocatable, dimension(:,:) :: edge_dist

    ! Set area cells array
    real(kind=8),allocatable, dimension(:) :: area

    integer :: bds
    real(kind=8) :: eps, fill_TH, marineslope, meanarea

    ! Sediment distribution parameters
    integer, parameter :: diff_nb = 2
    integer, parameter :: max_it_cyc = 10000

    real(kind=8), parameter :: diff_res = 1.e-4
    real(kind=8), parameter :: toplimit = 1.e10

contains

    subroutine defineparameters

      if(allocated(neighbours)) deallocate(neighbours)
      if(allocated(edge_dist)) deallocate(edge_dist)
      if(allocated(area)) deallocate(area)
      if(allocated(maxngb)) deallocate(maxngb)
      allocate(neighbours(dnodes,20))
      allocate(edge_dist(dnodes,20))
      allocate(area(dnodes))
      allocate(maxngb(dnodes))

      if(allocated(data1)) deallocate(data1)
      if(allocated(data2)) deallocate(data2)
      allocate(data1(block_size))
      allocate(data2(block_size))

      return

    end subroutine defineparameters

    subroutine find_maximum_elevation(difo, elev, sealevel, border, cdif, depo)

      integer :: i, p
      integer,dimension(dnodes),intent(in) :: border

      real(kind=8) :: topnew, topmax, sealevel
      real(kind=8),dimension(dnodes),intent(inout) :: cdif, depo
      real(kind=8),dimension(dnodes),intent(inout) :: difo, elev
      real(kind=8),dimension(dnodes) :: toph

      toph = elev
      do i = 1, dnodes
        if(difo(i)>0..and.area(i)>0..and.border(i)==1)then
          topmax = 1.0e6
          ! Get the maximum topographic elevation for each node
          loop: do p = 1, 20
            if(neighbours(i,p)<0) exit loop
            if(elev(i)<sealevel)then
             topnew = elev(neighbours(i,p)+1) + marineslope*edge_dist(neighbours(i,p)+1,p)
             if(topnew>sealevel)then
               topnew = sealevel + 0.0001*(edge_dist(neighbours(i,p)+1,p)+ &
                        (elev(i)-sealevel)/marineslope)
             endif
            else
             topnew = elev(neighbours(i,p)+1) + 0.0001*edge_dist(neighbours(i,p)+1,p)
            endif
            topmax = min(topmax,topnew)
          enddo loop

          ! Fraction of buffer that will be deposited at current node
          if(topmax>elev(i))then
            topnew = elev(i)+difo(i)/area(i)
            if(topnew<=topmax)then
              cdif(i) = 1.
              toph(i) = topnew
            else
              cdif(i) = (topmax-elev(i))/(topnew-elev(i))
              toph(i) = topmax
            endif
          else
            cdif(i) = 0.
            toph(i) = elev(i)
          endif
          depo(i) = depo(i)+difo(i)*cdif(i)/area(i)
        else
          difo(i)=0.
        endif
      enddo

      ! Update diffusion grid top elevation
      elev = toph

      return

    end subroutine find_maximum_elevation

    subroutine get_available_volume_remaining(k, elev, volsum, vol, down, ndown)

      integer,intent(in) :: k
      integer :: p
      integer, intent(out) :: ndown
      integer,dimension(20),intent(out) :: down

      real(kind=8) :: dz, slp
      real(kind=8),intent(out) :: volsum
      real(kind=8),dimension(20),intent(out) :: vol

      real(kind=8),dimension(dnodes),intent(in) :: elev

      volsum = 0.
      ndown = 0
      down = 0
      vol = 0.

      ! Loop over neighboring cell
      loop: do p = 1, 20
        if(neighbours(k,p)<0) exit loop
        dz = elev(k)-elev(neighbours(k,p)+1)
        if(dz>0.)then
          slp = dz/edge_dist(neighbours(k,p)+1,p)
          if(slp>=0.00001)then
            ndown = ndown+1
            down(p) = 1
            vol(p) = dz*area(neighbours(k,p)+1)
            volsum = volsum + vol(p)
          endif
        endif
      enddo loop

      return

    end subroutine get_available_volume_remaining

end module pdclass
