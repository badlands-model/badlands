!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Loop over node coordinates and find if they belong to local partition.

module part

  implicit none

contains

    subroutine overlap(pyX,pyY,pyXst,pyYst,pyXed,pyYed,pyPart,pyNodes)
     
        integer :: q
        integer :: pyNodes
        real(kind=8),intent(in) :: pyXst
        real(kind=8),intent(in) :: pyXed
        real(kind=8),intent(in) :: pyYst
        real(kind=8),intent(in) :: pyYed
        real(kind=8),intent(in) :: pyX(pyNodes)
        real(kind=8),intent(in) :: pyY(pyNodes)
        integer,intent(out) :: pyPart(pyNodes)
        
        pyPart = -1
        do q = 1,pyNodes
            if( pyX(q) >= pyXst .and. pyX(q) <= pyXed )then
                if( pyY(q) >= pyYst .and. pyY(q) <= pyYed )then
                    pyPart(q) = q-1
                endif
            endif
        enddo
        
        return

    end subroutine overlap

end module part