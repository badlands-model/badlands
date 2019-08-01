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
!       Filename:  Sorting.f90
!
!    Description:  Compute morphometrics based on hydrological features
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
! NOTE:
! This module arranges array elements from smallest to largest.
! REF:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! =====================================================================================
module sorting


  use parameters

  implicit none

  save

  private

  public quick_sort

contains

  ! =====================================================================================
  recursive subroutine quick_sort(list,order)

    real(kind=8),dimension(:),intent(inout)::list
    integer,dimension(:),intent(inout)::order

    call quick_sort_rec(1,size(list))

    contains

      ! =====================================================================================
      recursive subroutine quick_sort_rec(left_end,right_end)

        integer,intent(in)::left_end,right_end

        ! Local variables
        integer::i,j,itemp
        real(kind=8)::reference,temp
        integer,parameter::max_simple_sort_size=6

        if(right_end<left_end+max_simple_sort_size)then
          ! Use interchange sort for small lists
          call interchange_sort(left_end,right_end)
        else
          ! Use partition ("quick") sort
          reference=list((left_end+right_end)/2)
          i=left_end-1;j=right_end+1

          do
            ! Scan list from left end until element >= reference is found
            do
              i=i+1
              if(list(i)>=reference) exit
            enddo
            ! Scan list from right end until element <= reference is found
            do
              j=j-1
              if(list(j)<=reference) exit
            enddo

            if(i<j)then
              ! Swap two out-of-order elements
              temp=list(i);list(i)=list(j);list(j)=temp
              itemp=order(i);order(i)=order(j);order(j)=itemp
            elseif(i==j)then
              i=i+1
              exit
            else
              exit
            endif
          enddo

          if(left_end<j) call quick_sort_rec(left_end,j)
          if(i<right_end) call quick_sort_rec(i,right_end)
        endif

      end subroutine quick_sort_rec
      ! =====================================================================================

      subroutine interchange_sort(left_end,right_end)

        integer,intent(in)::left_end,right_end

        !  Local variables
        integer::i,j,itemp
        real(kind=8)::temp

        do i=left_end,right_end-1
          do j=i+1,right_end
            if(list(i)>list(j))then
              temp=list(i);list(i)=list(j);list(j)=temp
              itemp=order(i);order(i)=order(j);order(j)=itemp
            endif
          enddo
        enddo

      end subroutine interchange_sort
      ! =====================================================================================
  end subroutine quick_sort
  ! ============================================================================

end module sorting
