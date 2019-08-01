! =====================================================================================
! BALAD (BAsin and LAndscape Dynamics)
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
!       Filename:  UnderworldCpl.f90
!
!    Description:  Underworld coupling strategy
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module underworld

  use parallel
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none

contains

  ! =====================================================================================
  subroutine SurfaceVTK

    integer::iunit,ios,k,p,n,m,uorb

    uorb=0
    ! Create the top surface
    iunit=81

    if(pet_id==0)then
      open(iunit,file=fudw,status="replace",action="write",iostat=ios)
      rewind(iunit)

      ! Create header
      write(iunit,'(a26)')'# vtk DataFile Version 2.0'
      write(iunit,'(a24)')'Balad to Underworld Data'
      write(iunit,'(a5)')'ASCII'
      write(iunit,'(a25)')'DATASET STRUCTURED_POINTS'
      write(iunit,'(a10,1x,i10,1x,i10,1x,i2)')'DIMENSIONS',nx,ny,1
      write(iunit,'(a6,1x,f12.3,1x,f12.3,1x,f12.3)')'ORIGIN',minx,miny,0.0
      write(iunit,'(a7,1x,f12.3,1x,f12.3,1x,f12.3)')'SPACING',dx,dx,1.0
      write(iunit,'(a10,1x,i10)')'POINT_DATA',nx*ny
      write(iunit,'(a26)')'SCALARS elevation double 1'
      write(iunit,'(a20)')'LOOKUP_TABLE default'

      m=0
      do k=1,ny+2
        do n=1,nx+2
          m=m+1
          if(n>1.and.n<nx+1.and.k>1.and.k<ny+2)then
            write(iunit,'(f12.3,1x)',advance='no')rcoordZ(m)
          elseif(n==nx+1.and.k>1.and.k<ny+2)then
            write(iunit,'(f12.3,1x)')rcoordZ(m)
          endif
        enddo
      enddo

      write(iunit,'(a21)')'SCALARS fieldID int 1'
      write(iunit,'(a20)')'LOOKUP_TABLE default'

      p=0
      do k=1,ny+2
        do n=1,nx+2
          if(n>1.and.n<nx+1.and.k>1.and.k<ny+2)then
            p=p+1
            write(iunit,'(i6,1x)',advance='no')p
          elseif(n==nx+1.and.k>1.and.k<ny+2)then
            p=p+1
            write(iunit,'(i6)')p
          endif
        enddo
      enddo
      close(iunit)

      ! Update/Create the maestro file
      open(iunit,file=maestro,status="replace",action="write",iostat=ios)
      rewind(iunit)

      if(uorb==0)then
         write(iunit,'(a1)')'U'
         write(iunit,'(a1)')' '
      else
         write(iunit,'(a1)')'L'
         write(iunit,'(a1)')' '
      endif
      close(iunit)
    endif

    return

  end subroutine SurfaceVTK
  ! =====================================================================================
  subroutine WaitStepCompletion

    integer::iu,ios
    character(len=1)::charac

    charac='U'

    if(pet_id==0)then
       do while(charac/='L')
          iu=79
          ! Read the maestro file
          open(iu,file=maestro,status="old",action="read",iostat=ios)
          rewind(iu)
          if(ios==0)then
            read(iu,'(a1)',iostat=ios) charac
            if(ios/=0) charac='U'
          endif
          close(iu)
          call Sleep(1)
       enddo
    endif

    call mpi_barrier(badlands_world,rc)

    return

  end subroutine WaitStepCompletion
  ! =====================================================================================

end module underworld
! =====================================================================================
