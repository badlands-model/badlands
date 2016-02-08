!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Finite Volume discretisation implementation class for unstructured grid.

Module m_mrgrnk
Integer, Parameter :: kdp = selected_real_kind(15)
public :: mrgrnk
private :: kdp
private :: R_mrgrnk, I_mrgrnk, D_mrgrnk
interface mrgrnk
  module procedure D_mrgrnk, R_mrgrnk, I_mrgrnk
end interface mrgrnk 
contains

Subroutine D_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real (kind=kdp) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine D_mrgrnk

Subroutine R_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine R_mrgrnk
Subroutine I_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Integer :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine I_mrgrnk
end module m_mrgrnk



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



module fvclass

  use orderpack

  implicit none

  real(kind=8)::dminX,dminY,dmaxX,dmaxY

  ! Copy from Python
  integer :: onodes
  integer :: gnodes
  integer :: dnodes
  integer :: vnodes
  integer :: delems
  integer :: dedges
  integer :: vedges
  integer,dimension(:),allocatable :: oIDs,gIDs
  real(kind=8),dimension(:),allocatable :: tX
  real(kind=8),dimension(:),allocatable :: tY
  integer,dimension(:,:),allocatable :: tElmt
  integer,dimension(:,:),allocatable :: tEdge
  real(kind=8),dimension(:),allocatable :: vX
  real(kind=8),dimension(:),allocatable :: vY
  integer,dimension(:,:),allocatable :: vEdge
  integer,dimension(:),allocatable :: partIDs
  integer,dimension(:),allocatable :: ghosts

  ! Allocatated locally
  integer,dimension(:),allocatable::mask

  ! Voronoi Cell Declaration Type
  type voronoi
     integer::vertexNb
     integer::border
     integer::btype
     integer::bpoint
     integer,dimension(20)::vertexID
     real(kind=8)::area
  end type voronoi
  type(voronoi),dimension(:),allocatable::vor

  ! Conforming Delaunay Triangulation Declaration Type
  type delaunay
     integer::ngbNb
     integer,dimension(20)::ngbID
     real(kind=8),dimension(20)::voronoi_edge
     real(kind=8),dimension(20)::distance
  end type delaunay
  type(delaunay),dimension(:),allocatable::tri

contains

  ! =====================================================================================
  ! =====================================================================================
  ! Interactive module functions
  ! =====================================================================================
  ! =====================================================================================

  subroutine get_data

    integer :: n
    !real(kind=8),dimension(:,:),allocatable::map
    real(kind=8),dimension(2,vnodes)::map

    ! Delaunay data
    if(allocated(tri)) deallocate(tri)
    allocate(tri(dnodes))

    ! Voronoi data
    !if(allocated(map)) deallocate(map)
    if(allocated(mask)) deallocate(mask)
    !allocate(map(2,vnodes))
    allocate(mask(vnodes))
    dminX = 1.e8
    dminY = 1.e8
    dmaxX = -1.e8
    dmaxY = -1.e8
    do n=1,dnodes
        dminX = min(dminX,tX(n))
        dminY = min(dminY,tY(n))
        dmaxX = max(dmaxX,tX(n))
        dmaxY = max(dmaxY,tY(n))
    enddo
    
    ! Mask duplicate voronoi points
    mask = [(1+(n-1), n=1,vnodes)]
    map(1,:)=anint(vX(:)*100000.0)/100000.0
    map(2,:)=anint(vY(:)*100000.0)/100000.0

    call findmap(map,mask)

    do n=1,vedges
       if(vEdge(n,1)>0) vEdge(n,1)=mask(vEdge(n,1))
       if(vEdge(n,2)>0) vEdge(n,2)=mask(vEdge(n,2))
    enddo

    !deallocate(map)

    return

  end subroutine get_data
  
  subroutine build_export_arrays_o(pyVarea,pyNgbs,pyVlenght,pyDlenght,pymaxNgbh)

    integer :: cell,c,o,idd

    !integer,dimension(dnodes) :: pyDngbNn
    integer :: pymaxNgbh
    real(kind=8),dimension(dnodes) :: pyVarea
    integer,dimension(dnodes,20) :: pyNgbs
    real(kind=8),dimension(dnodes,20) :: pyVlenght
    real(kind=8),dimension(dnodes,20) :: pyDlenght
    
    pyVarea = 0.
    pyNgbs = -2
    !pyDngbNn = 0
    pyVlenght = 0.
    pyDlenght = 0.
    do c = 1, gnodes
      cell = gIDs(c)
      pyVarea(cell) = vor(cell)%area
      !pyDngbNn(cell) = tri(cell)%ngbNb
      pymaxNgbh = max(tri(cell)%ngbNb,pymaxNgbh)
      
      do o = 1,tri(cell)%ngbNb
          if(tri(cell)%ngbID(o)>0)then
              idd = oIDs(tri(cell)%ngbID(o))
              pyNgbs(cell,o) = idd-1
          endif
          pyVlenght(cell,o) = tri(cell)%voronoi_edge(o)
          pyDlenght(cell,o) = tri(cell)%distance(o)
      enddo
    enddo
    
  end subroutine build_export_arrays_o

  subroutine clean_data

    if(allocated(tri)) deallocate(tri)
    if(allocated(vor)) deallocate(vor)
    if(allocated(mask)) deallocate(mask)

    return

  end subroutine clean_data

  ! =====================================================================================
  ! =====================================================================================
  ! Delaunay Voronoi duality mesh
  ! =====================================================================================
  ! =====================================================================================

  subroutine delaunay_voronoi_duality(overOn)

    logical::rec,rec2,overOn

    integer::cell,n,id,id1,id2,id3,p,k,nvert,pt(2),c,k1,k2,pp
    integer,dimension(20)::vsort,sortedID

    !integer,dimension(:),allocatable::ed1,ed2,rk1,rk2
    integer,dimension(dedges)::ed1,ed2,rk1,rk2

    real(kind=8)::dist,area,s1,s2
    real(kind=8),dimension(2)::pt1,pt2,pt3,pt4
    real(kind=8),dimension(20)::vnx,vny

    ! Delaunay and voronoi parametrisation
    if(allocated(vor)) deallocate(vor)
    allocate(vor(dnodes))

    ! Rank the edges by point IDs
    !if(allocated(ed1)) deallocate(ed1,ed2)
    !if(allocated(rk1)) deallocate(rk1,rk2)
    !allocate(ed1(dedges),ed2(dedges))
    !allocate(rk1(dedges),rk2(dedges))
    ed1(:)=tEdge(:,1)
    ed2(:)=tEdge(:,2)
    call mrgrnk(ed1,rk1)
    call mrgrnk(ed2,rk2)
    id1=0
    id2=0
    k1=0
    k2=0
    
    do c = 1, gnodes
       cell = gIDs(c)
       id1=-1
       id2=-1
       lp_ed:do n=k1+1,dedges
         if(tEdge(rk1(n),1)==cell)then
           id1=n
           k1=n
           exit lp_ed
         endif
       enddo lp_ed
       
       lp_ed0:do n=k2+1,dedges
         if(tEdge(rk2(n),2)==cell)then
           id2=n
           k2=n
           exit lp_ed0
         endif
       enddo lp_ed0
       
       tri(cell)%ngbNb=0
       tri(cell)%ngbID=-1
       tri(cell)%voronoi_edge=0.0
       vor(cell)%vertexID=-1
       vor(cell)%vertexNb=0
       vor(cell)%border=0
       vor(cell)%area=0.0
       
       lp_ed1:do n=id1,dedges
          if(tEdge(rk1(n),1)==cell)then
             id1=n
             rec=.true.
             if(tEdge(rk1(n),2)>0)then
                ! Check that the points has not been recorded yet
                do id=1,tri(cell)%ngbNb
                   if(tri(cell)%ngbID(id)==tEdge(rk1(n),2)) rec=.false.
                enddo
                if(rec)then
                   id=tri(cell)%ngbNb+1
                   tri(cell)%ngbNb=id
                   tri(cell)%ngbID(id)=tEdge(rk1(n),2)
                endif
             endif

             ! Get the Voronoi cell by duality
             rec=.true.
             rec2=.true.
             if(vor(cell)%vertexNb>=1)then
                do p=1,vor(cell)%vertexNb
                   if(vor(cell)%vertexID(p)==vEdge(rk1(n),1)) rec=.false.
                   if(vor(cell)%vertexID(p)==vEdge(rk1(n),2)) rec2=.false.
                enddo
                if(rec.and.vEdge(rk1(n),1)>0)then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk1(n),1)
                endif
                if(rec2.and.vEdge(rk1(n),2)>0.and.vEdge(rk1(n),2)/=vEdge(rk1(n),1))then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk1(n),2)
                endif
             else
                if(vEdge(rk1(n),1)>0)then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk1(n),1)
                endif
                if(vEdge(rk1(n),2)>0.and.vEdge(rk1(n),2)/=vEdge(rk1(n),1))then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk1(n),2)
                endif
             endif
             if(vEdge(rk1(n),1)<=0) vor(cell)%border=1
          else
             exit lp_ed1
          endif
       enddo lp_ed1

       lp_ed2:do n=id2,dedges
          if(tEdge(rk2(n),2)==cell)then
             id2=n
             rec=.true.
             if(tEdge(rk2(n),1)>0)then
                ! Check that the points has not been recorded yet
                do id=1,tri(cell)%ngbNb
                   if(tri(cell)%ngbID(id)==tEdge(rk2(n),1)) rec=.false.
                enddo
                if(rec)then
                   id=tri(cell)%ngbNb+1
                   tri(cell)%ngbNb=id
                   tri(cell)%ngbID(id)=tEdge(rk2(n),1)
                endif
             endif

             ! Get the Voronoi cell by duality
             rec=.true.
             rec2=.true.
             if(vor(cell)%vertexNb>=1)then
                do p=1,vor(cell)%vertexNb
                   if(vor(cell)%vertexID(p)==vEdge(rk2(n),1)) rec=.false.
                   if(vor(cell)%vertexID(p)==vEdge(rk2(n),2)) rec2=.false.
                enddo
                if(rec.and.vEdge(rk2(n),1)>0)then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk2(n),1)
                endif
                if(rec2.and.vEdge(rk2(n),2)>0.and.vEdge(rk2(n),2)/=vEdge(rk2(n),1))then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk2(n),2)
                endif
             else
                if(vEdge(rk2(n),1)>0)then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk2(n),1)
                endif
                if(vEdge(rk2(n),2)>0.and.vEdge(rk2(n),2)/=vEdge(rk2(n),1))then
                   id=vor(cell)%vertexNb+1
                   vor(cell)%vertexNb=id
                   vor(cell)%vertexID(id)=vEdge(rk2(n),2)
                endif
             endif
             if(vEdge(rk2(n),1)<=0) vor(cell)%border=1
          else
             exit lp_ed2
          endif
       enddo lp_ed2
      
       ! Define the voronoi polygon as a convex hull
       if(vor(cell)%border==0)then
          vsort=-1
          n=vor(cell)%vertexNb
          do p=1,n
             id=vor(cell)%vertexID(p)
             vnx(p)=vX(id)
             vny(p)=vY(id)
          enddo

          ! Perform vertex sorting
          call envelope(vnx(1:n),vny(1:n),n,vsort,nvert)

          ! Reorganise voronoi cell vertex
          if(nvert/=vor(cell)%vertexNb)then
             if(tX(cell)>dminX.and.tX(cell)<dmaxX.and. &
                 tY(cell)>dminY.and.tY(cell)<dmaxY)then
                 print*,'Failed during creation of voronoi convex hull',cell
                 exit
             else
                 vor(cell)%border=1
             endif
          endif
          if(vor(cell)%border==0)then
              do p=1,nvert
                 sortedID(p)=vor(cell)%vertexID(vsort(p))
              enddo
              ! Update sorted vertex
              do p=1,nvert
                 vor(cell)%vertexID(p)=sortedID(p)
              enddo
              ! Get voronoi cell parameters
              do p=1,nvert
                 id=sortedID(p)
                 id3=sortedID(p+1)
                 if(p==nvert) id3=sortedID(1)
                 ! Area
                 area=(vX(id)-tX(cell))*(vY(id3)-tY(cell)) &
                      -(vX(id3)-tX(cell))*(vY(id)-tY(cell))
                 vor(cell)%area=vor(cell)%area+0.5*abs(area)
              enddo
          endif
       endif   
    enddo

    ! Allocate voronoi edge length to delaunay neighborhood
    pp = 1
    do c = 1, gnodes
       cell = gIDs(c)
       if(vor(cell)%border==0)then
          pt1(1)=tX(cell)
          pt1(2)=tY(cell)
          do n=1,tri(cell)%ngbNb
             id=tri(cell)%ngbID(n)
             if(.not. overOn)then
                 if( partIDs(cell)/=partIDs(id))then 
                     ghosts(pp)=id-1
                     pp=pp+1
                 endif
             endif
             ! Find voronoi edge corresponding to this neighbour
             pt2(1)=tX(id)
             pt2(2)=tY(id)
             s1=0._8
             pt = -1
             if(pt1(1)/=pt2(1)) s1=abs((pt2(2)-pt1(2))/(pt2(1)-pt1(1)))
             llpp:do p=1,vor(cell)%vertexNb
                k=p+1
                if(p==vor(cell)%vertexNb) k=1
                pt3(1)=vX(vor(cell)%vertexID(p))
                pt3(2)=vY(vor(cell)%vertexID(p))
                pt4(1)=vX(vor(cell)%vertexID(k))
                pt4(2)=vY(vor(cell)%vertexID(k))
                s2=0._8
                if(pt4(1)/=pt3(1)) s2=abs((pt3(2)-pt4(2))/(pt3(1)-pt4(1)))
                if(abs(s1*s2-1.0_8)<0.001_8)then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=vor(cell)%vertexID(k)
                        pt(2)=vor(cell)%vertexID(p)
                        exit llpp
                    endif
                else if(pt1(1)==pt2(1).and.pt4(2)==pt3(2))then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=vor(cell)%vertexID(k)
                        pt(2)=vor(cell)%vertexID(p)
                        exit llpp
                    endif
                else if(pt1(2)==pt2(2).and.pt4(1)==pt3(1))then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=vor(cell)%vertexID(k)
                        pt(2)=vor(cell)%vertexID(p)
                        exit llpp
                    endif
                endif
             enddo llpp
             dist=0.
             if(pt(1)>0.and.pt(2)>0) &
                 dist=sqrt((vX(pt(1))-vX(pt(2)))**2.0+(vY(pt(2))-vY(pt(1)))**2.0)
             tri(cell)%voronoi_edge(n)=dist
             tri(cell)%distance(n)=sqrt((tX(cell)-tX(id))**2.0+( &
                  tY(cell)-tY(id))**2.0)
          enddo
       endif
    enddo
    
    !if(allocated(ed1)) deallocate(ed1,ed2)
    !if(allocated(rk1)) deallocate(rk1,rk2)

    return

  end subroutine delaunay_voronoi_duality

  ! =====================================================================================
  ! =====================================================================================
  ! Usefull geometrical functions
  ! =====================================================================================
  ! =====================================================================================

  subroutine envelope(x,y,n,vertex,nvert)

    ! Find the vertices (in clockwise order) of a polygon enclosing
    ! the points (x(i), y(i), i=1, ..., n.
    ! On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
    ! iwk() is an integer work array which must have dimension at least n
    ! in the calling program.

    integer::n,vertex(n),nvert,iwk(n)
    integer::next(20),i,i1,i2,j,jp1,jp2,i2save,i3,i2next

    real(kind=8)::x(n),y(n),xmax,xmin,ymax,ymin,dist,dmax,dmin,x1,y1
    real(kind=8)::dx,dy,x2,y2,dx1,dx2,dmax1,dmax2,dy1,dy2,temp,zero

    zero=0.0_8
    x = x
    y = y
    if(n<2) return

    ! Choose the points with smallest & largest x- values as the
    ! first two vertices of the polygon.

    if(x(1)>x(n))then
       vertex(1)=n
       vertex(2)=1
       xmin=x(n)
       xmax=x(1)
    else
       vertex(1)=1
       vertex(2)=n
       xmin=x(1)
       xmax=x(n)
    endif

    next = -1
    do i=2,n-1
       temp=x(i)
       if(temp<xmin)then
          vertex(1)=i
          xmin=temp
       elseif(temp>xmax)then
          vertex(2)=i
          xmax=temp
       endif
    enddo

    ! Special case, xmax = xmin
    if(xmax==xmin)then
       if(y(1)>y(n))then
          vertex(1)=n
          vertex(2)=1
          ymin=y(n)
          ymax=y(1)
       else
          vertex(1)=1
          vertex(2)=n
          ymin=y(1)
          ymax=y(n)
       endif

       do i=2,n-1
          temp=y(i)
          if(temp<ymin)then
             vertex(1)=i
             ymin=temp
          elseif(temp>ymax)then
             vertex(2)=i
             ymax=temp
          endif
       enddo

       nvert=2
       if(ymax==ymin) nvert=1
       return
    endif

    ! Set up two initial lists of points; those points above & those below the
    ! line joining the first two vertices. next(i) will hold the pointer to the
    ! point furthest from the line joining vertex(i) to vertex(i+1) on the left
    ! hand side.
    i1=vertex(1)
    i2=vertex(2)
    iwk(i1)=-1
    iwk(i2)=-1
    dx=xmax-xmin
    y1=y(i1)
    dy=y(i2)-y1
    dmax=zero
    dmin=zero
    next(1)=-1
    next(2)=-1

    do i=1,n
       if(i==vertex(1) .or. i==vertex(2)) cycle
       dist=(y(i)-y1)*dx-(x(i)-xmin)*dy
       if(dist>zero)then
          iwk(i1)=i
          i1=i
          if(dist>dmax)then
             next(1)=i
             dmax=dist
          endif
       elseif(dist<zero)then
          iwk(i2)=i
          i2=i
          if(dist<dmin)then
             next(2)=i
             dmin=dist
          endif
       endif
    enddo

    ! Ends of lists are indicated by pointers to -ve positions.
    iwk(i1)=-1
    iwk(i2)=-1
    nvert=2

    j=1

    ! Start of main process.
    ! Introduce new vertex between vertices j & j+1, if one has been found.
    ! Otherwise increase j. Exit if no more vertices.
40  if(next(j)<0)then
       if(j==nvert) return
       j=j+1
       goto 40
    endif

    jp1=j+1
    do i=nvert,jp1,-1
       vertex(i+1)=vertex(i)
       next(i+1)=next(i)
    enddo
    jp2=jp1+1
    nvert=nvert+1
    if(jp2>nvert) jp2=1
    i1=vertex(j)
    i2=next(j)
    i3=vertex(jp2)
    vertex(jp1)=i2

    ! Process the list of points associated with vertex j.   New list at vertex j
    ! consists of those points to the left of the line joining it to the new
    ! vertex (j+1).   Similarly for the list at the new vertex.
    ! Points on or to the right of these lines are dropped.
    x1=x(i1)
    x2=x(i2)
    y1=y(i1)
    y2=y(i2)
    dx1=x2-x1
    dx2=x(i3)-x2
    dy1=y2-y1
    dy2=y(i3)-y2
    DMAX1=zero
    dmax2=zero
    next(j)=-1
    next(jp1)=-1
    i2save=i2
    i2next=iwk(i2)
    i=iwk(i1)
    iwk(i1)=-1
    iwk(i2)=-1

60  if(i/=i2save)then
       dist=(y(i)-y1)*dx1-(x(i)-x1)*dy1
       if(dist>zero)then
          iwk(i1)=i
          i1=i
          if(dist>DMAX1)then
             next(j)=i
             DMAX1=dist
          endif
       else
          dist=(y(i)-y2)*dx2-(x(i)-x2)*dy2
          if(dist>zero)then
             iwk(i2)=i
             i2=i
             if(dist>dmax2)then
                next(jp1)=i
                dmax2=dist
             endif
          endif
       endif
       i=iwk(i)
    else
       i=i2next
    endif

    ! Get next point from old list at vertex j.
    if(i>0) goto 60

    ! End lists with -ve values.
    iwk(i1)=-1
    iwk(i2)=-1

    goto 40

  end subroutine envelope
 
  function find_intersection_between_segments(pt1,pt2,pt3,pt4) result(solution)

      logical :: solution
      real(kind=8),dimension(2)::pt1,pt2,pt3,pt4
      real(kind=8),dimension(2)::p13,p43,p21 !,intersect
      real(kind=8)::num_a,num_b,denom,mua,mub
      !integer::sol

      solution=.false.
      !intersect(:)=-1.e7_8
      p21(1)=pt2(1)-pt1(1)
      p21(2)=pt2(2)-pt1(2)
      p43(1)=pt4(1)-pt3(1)
      p43(2)=pt4(2)-pt3(2)
      p13(1)=pt1(1)-pt3(1)
      p13(2)=pt1(2)-pt3(2)
      denom=p43(2)*p21(1)-p43(1)*p21(2)
      num_a=p43(1)*p13(2)-p43(2)*p13(1)
      num_b=p21(1)*p13(2)-p21(2)*p13(1)
      if(denom==0.0_8)then
         return
      else
         mua=num_a/denom
         mub=num_b/denom
         ! Intersecting
         if(mua>=0.0_8 .and. mua<=1.0_8)then
            solution = .true.
            return
         else
            return
         endif   
      endif

      return

   end function find_intersection_between_segments
   
end module fvclass
