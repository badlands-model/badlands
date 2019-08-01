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
!       Filename:  Orderpack.f90
!
!    Description:  Find duplicate records in fortran
!
!        Version:  1.0
!        Created:  29/04/15 08:18:31
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module orderpack

  use m_mrgrnk, only:mrgrnk

  implicit none

contains
  ! =====================================================================================
  subroutine findmap(stkprm,stkmap)

    ! Given 2-d real array stkprm, find a mapping described below:
    !
    ! (identical records are assigned with same index)
    !   stkmap(i) == stkmap(j)  if stkprm(:,i) == stkprm(:,j)
    ! (order conserved)
    !   if i < j and stkmap(i) /= stkmap(j), then stkmap(i) < stkmap(j)
    ! (new index are contiguous)
    !   set(stkmap) == {1,2,..,maxval(stkmap)}
    !
    real(kind=8),dimension(:,:),intent(in) :: stkprm
!     integer,dimension(:) :: stkmap
    integer,dimension(:), intent(out) :: stkmap
    integer, dimension(size(stkprm,2)) :: irngt
    integer, dimension(size(stkprm,2)) :: iwork
    integer ::  nrec, i, j

    nrec = size(stkprm,2)

    ! Find rank of each record, duplicate records kept
    call ar_mrgrnk(stkprm, irngt)

    ! construct iwork array, which has index of original array where the
    ! record are identical, and the index is younguest
    i = 1
    do while(i<=nrec)
      do j=i+1,nrec
        if (any(stkprm(:,irngt(i))/=stkprm(:,irngt(j)))) exit
      enddo
      iwork(irngt(i:j-1)) = minval(irngt(i:j-1))
      i = j
    enddo

    stkmap=iwork

    ! now construct the map, where stkmap(i) shows index of new array
    ! with duplicated record eliminated, original order kept
!     j = 0
!     do i=1,nrec
!       if (i==iwork(i)) then
!         j = j+1
!         stkmap(i) = j
!       else
!         stkmap(i) = stkmap(iwork(i))
!         print*,'dde',i,iwork(i),stkmap(i),stkmap(iwork(i))
!       endif
!     enddo

  end subroutine findmap
  ! =====================================================================================
  recursive subroutine ar_mrgrnk(xdont, irngt)

    ! behaves like mrgrnk of ORDERPACK, except that array is 2-d
    ! each row are ranked by first field, then second and so on
    real(kind=8), dimension(:,:), intent(in) :: xdont
    integer, dimension(:), intent(out), target :: irngt
    integer, dimension(size(xdont,2)) :: iwork

    integer :: nfld,nrec
    integer :: i, j
    integer, dimension(:), pointer :: ipt

    nfld=size(xdont,1)
    nrec=size(xdont,2)

    ! rank by the first field
    call mrgrnk(xdont(1,:), irngt)

    ! if there's only one field, it's done
    if (nfld==1) return

    ! examine the rank to see if multiple record has identical
    ! values for the first field
    i = 1
    do while(i<=nrec)
      do j=i+1,nrec
        if (xdont(1,irngt(i))/=xdont(1,irngt(j))) exit
      enddo
      ! if one-to-one, do nothing
      if (j-1>i) then
      ! if many-to-one,
        ! gather those many, and rank them
        call ar_mrgrnk(xdont(2:,irngt(i:j-1)),iwork)
        ! rearrange my rank based on those fields to the right
        ipt => irngt(i:j-1)
        ipt = ipt(iwork(1:j-i))
      endif
      i = j
    enddo
    if(associated(ipt)) nullify(ipt)

  end subroutine ar_mrgrnk
  ! =====================================================================================
end module orderpack
! =====================================================================================
