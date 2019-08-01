!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Finite Volume discretisation implementation class for unstructured grid.

module m_mrgrnk

  integer, Parameter :: kdp = selected_real_kind(15)
  public :: mrgrnk
  private :: kdp
  private :: R_mrgrnk, I_mrgrnk, D_mrgrnk

  interface mrgrnk
    module procedure D_mrgrnk, R_mrgrnk, I_mrgrnk
  end interface mrgrnk

contains

  ! =====================================================================================
  subroutine D_mrgrnk (XdoNT, IRNGT)

      real (kind=kdp), dimension (:), intent (in) :: XdoNT
      integer, dimension (:), intent (out) :: IRNGT
      real (kind=kdp) :: XVALA, XVALB
      integer, dimension (size(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IinD, IWRKD, IWRK, IWRKF, JinDA, IinDA, IinDB

      NVAL = min (size(XdoNT), size(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select

      do IinD = 2, NVAL, 2
         if (XdoNT(IinD-1) <= XdoNT(IinD)) then
            IRNGT (IinD-1) = IinD - 1
            IRNGT (IinD) = IinD
         else
            IRNGT (IinD-1) = IinD
            IRNGT (IinD) = IinD - 1
         endif
      enddo
      if (modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      endif

      LMTNA = 2
      LMTNC = 4

      do
         if (NVAL <= 2) exit
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
               if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) exit
               if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               endif
               exit
            endif
            if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) Cycle
            if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               endif
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG1) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                     IRNGT (IWRKD+3) = IRNG2
                  else
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  endif
               else
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               endif
            endif
         enddo
         LMTNA = 4
         exit
      enddo
      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JinDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JinDA >= NVAL) exit
               IWRKF = NVAL
            endif
            IinDA = 1
            IinDB = JinDA + 1
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JinDA)
            XVALA = XdoNT (JWRKT(IinDA))
            XVALB = XdoNT (IRNGT(IinDB))
            do
               IWRK = IWRK + 1
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IinDB)
                  IinDB = IinDB + 1
                  if (IinDB > IWRKF) then
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IinDA:LMTNA)
                     exit
                  endif
                  XVALB = XdoNT (IRNGT(IinDB))
               else
                  IRNGT (IWRK) = JWRKT (IinDA)
                  IinDA = IinDA + 1
                  if (IinDA > LMTNA) exit
                  XVALA = XdoNT (JWRKT(IinDA))
               endif
            enddo
         enddo
         LMTNA = 2 * LMTNA
      enddo

      return

  end subroutine D_mrgrnk
  ! =====================================================================================
  subroutine R_mrgrnk (XdoNT, IRNGT)

      real, dimension (:), intent (in) :: XdoNT
      integer, dimension (:), intent (out) :: IRNGT
      real :: XVALA, XVALB
      integer, dimension (size(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IinD, IWRKD, IWRK, IWRKF, JinDA, IinDA, IinDB

      NVAL = min (size(XdoNT), size(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select
      do IinD = 2, NVAL, 2
         if (XdoNT(IinD-1) <= XdoNT(IinD)) then
            IRNGT (IinD-1) = IinD - 1
            IRNGT (IinD) = IinD
         else
            IRNGT (IinD-1) = IinD
            IRNGT (IinD) = IinD - 1
         endif
      enddo
      if (modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      endif
      LMTNA = 2
      LMTNC = 4
      do
         if (NVAL <= 2) exit
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
               if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) exit
               if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               endif
               exit
            endif
            if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) Cycle
            if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               endif
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG1) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                     IRNGT (IWRKD+3) = IRNG2
                  else
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  endif
               else
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               endif
            endif
         enddo
         LMTNA = 4
         exit
      enddo

      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JinDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JinDA >= NVAL) exit
               IWRKF = NVAL
            endif
            IinDA = 1
            IinDB = JinDA + 1
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JinDA)
            XVALA = XdoNT (JWRKT(IinDA))
            XVALB = XdoNT (IRNGT(IinDB))
            do
               IWRK = IWRK + 1
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IinDB)
                  IinDB = IinDB + 1
                  if (IinDB > IWRKF) then
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IinDA:LMTNA)
                     exit
                  endif
                  XVALB = XdoNT (IRNGT(IinDB))
               else
                  IRNGT (IWRK) = JWRKT (IinDA)
                  IinDA = IinDA + 1
                  if (IinDA > LMTNA) exit
                  XVALA = XdoNT (JWRKT(IinDA))
               endif
            enddo
         enddo
         LMTNA = 2 * LMTNA
      enddo

      return
  end subroutine R_mrgrnk
  ! =====================================================================================
  subroutine I_mrgrnk (XdoNT, IRNGT)

      integer, dimension (:), intent (in)  :: XdoNT
      integer, dimension (:), intent (out) :: IRNGT
      integer :: XVALA, XVALB
      integer, dimension (size(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IinD, IWRKD, IWRK, IWRKF, JinDA, IinDA, IinDB

      NVAL = min (size(XdoNT), size(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select
      do IinD = 2, NVAL, 2
         if (XdoNT(IinD-1) <= XdoNT(IinD)) then
            IRNGT (IinD-1) = IinD - 1
            IRNGT (IinD) = IinD
         else
            IRNGT (IinD-1) = IinD
            IRNGT (IinD) = IinD - 1
         endif
      enddo
      if (modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      endif
      LMTNA = 2
      LMTNC = 4
      do
         if (NVAL <= 2) exit
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
               if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) exit
               if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               endif
               exit
            endif
            if (XdoNT(IRNGT(IWRKD+2)) <= XdoNT(IRNGT(IWRKD+3))) Cycle
            if (XdoNT(IRNGT(IWRKD+1)) <= XdoNT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+3) = IRNG2
               else
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               endif
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XdoNT(IRNG1) <= XdoNT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XdoNT(IRNG2) <= XdoNT(IRNGT(IWRKD+4))) then
                     IRNGT (IWRKD+3) = IRNG2
                  else
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  endif
               else
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               endif
            endif
         enddo
         LMTNA = 4
         exit
      enddo

      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JinDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JinDA >= NVAL) exit
               IWRKF = NVAL
            endif
            IinDA = 1
            IinDB = JinDA + 1
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JinDA)
            XVALA = XdoNT (JWRKT(IinDA))
            XVALB = XdoNT (IRNGT(IinDB))
            do
               IWRK = IWRK + 1
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IinDB)
                  IinDB = IinDB + 1
                  if (IinDB > IWRKF) then
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IinDA:LMTNA)
                     exit
                  endif
                  XVALB = XdoNT (IRNGT(IinDB))
               else
                  IRNGT (IWRK) = JWRKT (IinDA)
                  IinDA = IinDA + 1
                  if (IinDA > LMTNA) exit
                  XVALA = XdoNT (JWRKT(IinDA))
               endif
            enddo
         enddo
         LMTNA = 2 * LMTNA
      enddo

      return

  end subroutine I_mrgrnk

end module m_mrgrnk
! =====================================================================================
! =====================================================================================
module orderpack

  use m_mrgrnk, only:mrgrnk

  implicit none

contains

  ! =====================================================================================
  subroutine findmap(stkprm,stkmap)

    real(kind=8),dimension(:,:),intent(in) :: stkprm
    integer,dimension(:), intent(out) :: stkmap
    integer, dimension(size(stkprm,2)) :: irngt
    integer, dimension(size(stkprm,2)) :: iwork
    integer ::  nrec, i, j

    nrec = size(stkprm,2)
    call ar_mrgrnk(stkprm, irngt)
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

    real(kind=8), dimension(:,:), intent(in) :: xdont
    integer, dimension(:), intent(out), target :: irngt
    integer, dimension(size(xdont,2)) :: iwork
    integer :: nfld,nrec
    integer :: i, j
    integer, dimension(:), pointer :: ipt

    nfld=size(xdont,1)
    nrec=size(xdont,2)
    call mrgrnk(xdont(1,:), irngt)
    if (nfld==1) return
    i = 1
    do while(i<=nrec)
      do j=i+1,nrec
        if (xdont(1,irngt(i))/=xdont(1,irngt(j))) exit
      enddo
      if (j-1>i) then
        call ar_mrgrnk(xdont(2:,irngt(i:j-1)),iwork)
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
! =====================================================================================
module classfv

  use orderpack

  implicit none

  real(kind=8)::dminX,dminY,dmaxX,dmaxY

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

  integer,dimension(:),allocatable::mask

  integer,dimension(:),allocatable :: voronoi_vertexNb
  integer,dimension(:),allocatable :: voronoi_border
  integer,dimension(:,:),allocatable :: voronoi_vertexID
  real(kind=8),dimension(:),allocatable :: voronoi_area

  integer,dimension(:),allocatable :: tri_ngbNb
  integer,dimension(:,:),allocatable :: tri_ngbID
  integer,dimension(:,:),allocatable :: tri_voronoi_edge
  real(kind=8),dimension(:,:),allocatable :: tri_distance

  integer :: incisiontype
  integer :: bedslptype
  real(kind=8) :: spl_n
  real(kind=8) :: spl_m
  real(kind=8) :: sed_mt
  real(kind=8) :: sed_nt
  real(kind=8) :: sed_kt
  real(kind=8) :: width_kw
  real(kind=8) :: width_b

  integer,dimension(:),allocatable :: allocs
  integer,dimension(:),allocatable :: Donors
  integer,dimension(:),allocatable :: Delta
  integer,dimension(:),allocatable :: stackOrder

contains

  recursive function addtostack(base,donor,stackID) result(success)

      integer :: base,donor,stackID,n,success

      success = 1
      do n = Delta(donor),Delta(donor+1)-1
          if((allocs(Donors(n)) /= base))then
              donor = Donors(n)
              stackID = stackID+1
              stackOrder(stackID) = donor
              allocs(donor) = base
              success = addtostack(base,donor,stackID)
          endif
      enddo
      success = 0

  end function addtostack
  ! =====================================================================================
  subroutine get_data

    integer :: n
    real(kind=8),dimension(2,vnodes)::map

    if(allocated(tri_ngbNb)) deallocate(tri_ngbNb)
    if(allocated(tri_ngbID)) deallocate(tri_ngbID)
    if(allocated(tri_distance)) deallocate(tri_distance)
    if(allocated(tri_voronoi_edge)) deallocate(tri_voronoi_edge)
    allocate(tri_ngbNb(dnodes))
    allocate(tri_ngbID(dnodes,20))
    allocate(tri_distance(dnodes,20))
    allocate(tri_voronoi_edge(dnodes,20))
    if(allocated(mask)) deallocate(mask)
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
    mask = [(1+(n-1), n=1,vnodes)]
    map(1,:)=anint(vX(:)*100000.0)/100000.0
    map(2,:)=anint(vY(:)*100000.0)/100000.0
    call findmap(map,mask)
    do n=1,vedges
       if(vEdge(n,1)>0) vEdge(n,1)=mask(vEdge(n,1))
       if(vEdge(n,2)>0) vEdge(n,2)=mask(vEdge(n,2))
    enddo

    return

  end subroutine get_data
  ! =====================================================================================
  subroutine build_export_arrays_o(pyVarea,pyNgbs,pyVlenght,pyDlenght,pymaxNgbh)

    integer :: cell,c,o,idd,o2
    integer :: pymaxNgbh
    real(kind=8),dimension(dnodes) :: pyVarea
    integer,dimension(dnodes,20) :: pyNgbs
    real(kind=8),dimension(dnodes,20) :: pyVlenght
    real(kind=8),dimension(dnodes,20) :: pyDlenght

    pyVarea = 0.
    pyNgbs = -2
    pyVlenght = 0.
    pyDlenght = 0.
    do c = 1, gnodes
      cell = gIDs(c)
      pyVarea(cell) = voronoi_area(cell)
      pymaxNgbh = max(tri_ngbNb(cell),pymaxNgbh)
      do o = 1,tri_ngbNb(cell)
          if(tri_ngbID(cell,o)>0)then
              idd = oIDs(tri_ngbID(cell,o))
              pyNgbs(cell,o) = idd-1
          endif
          pyVlenght(cell,o) = tri_voronoi_edge(cell,o)
          pyDlenght(cell,o) = tri_distance(cell,o)
      enddo
    enddo

    do c = 1, gnodes
      cell = gIDs(c)
      do o = 1,tri_ngbNb(cell)
        if(tri_ngbID(cell,o)>0)then
            idd = pyNgbs(cell,o)+1
            do o2 = 1,tri_ngbNb(idd)
              if(pyNgbs(idd,o2) == cell-1)then
                if(pyDlenght(cell,o).ne.pyDlenght(idd,o2))then
                  if(pyDlenght(cell,o)==0.)pyDlenght(cell,o) = pyDlenght(idd,o2)
                  if(pyDlenght(idd,o2)==0.)pyDlenght(idd,o2) = pyDlenght(cell,o)
                endif
                if(pyVlenght(cell,o).ne.pyVlenght(idd,o2))then
                  pyVlenght(cell,o) = pyVlenght(idd,o2)
                endif
              endif
            enddo
        endif
      enddo
    enddo

    return

  end subroutine build_export_arrays_o
  ! =====================================================================================
  subroutine delaunay_voronoi_duality(overOn)

    logical::rec,rec2,overOn
    integer::cell,n,id,id1,id2,id3,p,k,nvert,pt(2),c,k1,k2,pp
    integer,dimension(20)::vsort,sortedID
    integer,dimension(dedges)::ed1,ed2,rk1,rk2
    real(kind=8)::dist,area,s1,s2
    real(kind=8),dimension(2)::pt1,pt2,pt3,pt4
    real(kind=8),dimension(20)::vnx,vny

    if(allocated(voronoi_area)) deallocate(voronoi_area)
    if(allocated(voronoi_border)) deallocate(voronoi_border)
    if(allocated(voronoi_vertexNb)) deallocate(voronoi_vertexNb)
    if(allocated(voronoi_vertexID)) deallocate(voronoi_vertexID)
    allocate(voronoi_area(dnodes))
    allocate(voronoi_border(dnodes))
    allocate(voronoi_vertexNb(dnodes))
    allocate(voronoi_vertexID(dnodes,20))
    ed1(:)=tEdge(:,1)
    ed2(:)=tEdge(:,2)
    call mrgrnk(ed1,rk1)
    call mrgrnk(ed2,rk2)
    id1=0
    id2=0
    k1=0
    k2=0

    tri_distance = 0.
    voronoi_area = 0.
    tri_voronoi_edge = 0.
    voronoi_vertexID = -1
    voronoi_vertexNb = 0
    voronoi_border = 0
    voronoi_area = 0.0
    tri_ngbNb = 0
    tri_ngbID = -1


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
       lp_ed0: do n=k2+1,dedges
         if(tEdge(rk2(n),2)==cell)then
           id2=n
           k2=n
           exit lp_ed0
         endif
       enddo lp_ed0

       lp_ed1: do n=id1,dedges
          if(tEdge(rk1(n),1)==cell)then
             id1=n
             rec=.true.
             if(tEdge(rk1(n),2)>0)then
                do id=1,tri_ngbNb(cell)
                   if(tri_ngbID(cell,id)==tEdge(rk1(n),2)) rec=.false.
                enddo
                if(rec)then
                   id=tri_ngbNb(cell)+1
                   tri_ngbNb(cell)=id
                   tri_ngbID(cell,id)=tEdge(rk1(n),2)
                endif
             endif
             rec=.true.
             rec2=.true.
             if(voronoi_vertexNb(cell)>=1)then
                do p=1,voronoi_vertexNb(cell)
                   if(voronoi_vertexID(cell,p)==vEdge(rk1(n),1)) rec=.false.
                   if(voronoi_vertexID(cell,p)==vEdge(rk1(n),2)) rec2=.false.
                enddo
                if(rec.and.vEdge(rk1(n),1)>0)then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk1(n),1)
                endif
                if(rec2.and.vEdge(rk1(n),2)>0.and.vEdge(rk1(n),2)/=vEdge(rk1(n),1))then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk1(n),2)
                endif
             else
                if(vEdge(rk1(n),1)>0)then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk1(n),1)
                endif
                if(vEdge(rk1(n),2)>0.and.vEdge(rk1(n),2)/=vEdge(rk1(n),1))then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk1(n),2)
                endif
             endif
             if(vEdge(rk1(n),1)<=0) voronoi_border(cell)=1
          else
             exit lp_ed1
          endif
       enddo lp_ed1

       lp_ed2: do n=id2,dedges
          if(tEdge(rk2(n),2)==cell)then
             id2=n
             rec=.true.
             if(tEdge(rk2(n),1)>0)then
                do id=1,tri_ngbNb(cell)
                   if(tri_ngbID(cell,id)==tEdge(rk2(n),1)) rec=.false.
                enddo
                if(rec)then
                   id=tri_ngbNb(cell)+1
                   tri_ngbNb(cell)=id
                   tri_ngbID(cell,id)=tEdge(rk2(n),1)
                endif
             endif
             rec=.true.
             rec2=.true.
             if(voronoi_vertexNb(cell)>=1)then
                do p=1,voronoi_vertexNb(cell)
                   if(voronoi_vertexID(cell,p)==vEdge(rk2(n),1)) rec=.false.
                   if(voronoi_vertexID(cell,p)==vEdge(rk2(n),2)) rec2=.false.
                enddo
                if(rec.and.vEdge(rk2(n),1)>0)then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk2(n),1)
                endif
                if(rec2.and.vEdge(rk2(n),2)>0.and.vEdge(rk2(n),2)/=vEdge(rk2(n),1))then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk2(n),2)
                endif
             else
                if(vEdge(rk2(n),1)>0)then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk2(n),1)
                endif
                if(vEdge(rk2(n),2)>0.and.vEdge(rk2(n),2)/=vEdge(rk2(n),1))then
                   id=voronoi_vertexNb(cell)+1
                   voronoi_vertexNb(cell)=id
                   voronoi_vertexID(cell,id)=vEdge(rk2(n),2)
                endif
             endif
             if(vEdge(rk2(n),1)<=0) voronoi_border(cell)=1
          else
             exit lp_ed2
          endif
       enddo lp_ed2

       if(voronoi_border(cell)==0)then
          vsort=-1
          n=voronoi_vertexNb(cell)
          do p=1,n
             id=voronoi_vertexID(cell,p)
             vnx(p)=vX(id)
             vny(p)=vY(id)
          enddo
          call envelope(vnx(1:n),vny(1:n),n,vsort,nvert)
          if(nvert/=voronoi_vertexNb(cell))then
             if(tX(cell)>dminX.and.tX(cell)<dmaxX.and. &
                 tY(cell)>dminY.and.tY(cell)<dmaxY)then
                 print*,'Failed during creation of voronoi convex hull',cell
                 exit
             else
                 voronoi_border(cell)=1
             endif
          endif
          if(voronoi_border(cell)==0)then
              do p=1,nvert
                 sortedID(p)=voronoi_vertexID(cell,vsort(p))
              enddo
              do p=1,nvert
                 voronoi_vertexID(cell,p)=sortedID(p)
              enddo
              do p=1,nvert
                 id=sortedID(p)
                 id3=sortedID(p+1)
                 if(p==nvert) id3=sortedID(1)
                 area=(vX(id)-tX(cell))*(vY(id3)-tY(cell)) &
                      -(vX(id3)-tX(cell))*(vY(id)-tY(cell))
                 voronoi_area(cell)=voronoi_area(cell)+0.5*abs(area)
              enddo
          endif
       endif
    enddo

    pp = 1
    do c = 1, gnodes
       cell = gIDs(c)
       if(voronoi_border(cell)==0)then
          pt1(1)=tX(cell)
          pt1(2)=tY(cell)
          do n=1,tri_ngbNb(cell)
             id=tri_ngbID(cell,n)
             if(.not. overOn)then
                 if( partIDs(cell)/=partIDs(id))then
                     ghosts(pp)=id-1
                     pp=pp+1
                 endif
             endif
             pt2(1)=tX(id)
             pt2(2)=tY(id)
             s1=0._8
             pt = -1
             if(pt1(1)/=pt2(1)) s1=abs((pt2(2)-pt1(2))/(pt2(1)-pt1(1)))
             llpp: do p=1,voronoi_vertexNb(cell)
                k=p+1
                if(p==voronoi_vertexNb(cell)) k=1
                pt3(1)=vX(voronoi_vertexID(cell,p))
                pt3(2)=vY(voronoi_vertexID(cell,p))
                pt4(1)=vX(voronoi_vertexID(cell,k))
                pt4(2)=vY(voronoi_vertexID(cell,k))
                s2=0._8
                if(pt4(1)/=pt3(1)) s2=abs((pt3(2)-pt4(2))/(pt3(1)-pt4(1)))
                if(abs(s1*s2-1.0_8)<0.001_8)then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=voronoi_vertexID(cell,k)
                        pt(2)=voronoi_vertexID(cell,p)
                        exit llpp
                    endif
                else if(pt1(1)==pt2(1).and.pt4(2)==pt3(2))then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=voronoi_vertexID(cell,k)
                        pt(2)=voronoi_vertexID(cell,p)
                        exit llpp
                    endif
                else if(pt1(2)==pt2(2).and.pt4(1)==pt3(1))then
                    if(find_intersection_between_segments(pt1,pt2,pt3,pt4))then
                        pt(1)=voronoi_vertexID(cell,k)
                        pt(2)=voronoi_vertexID(cell,p)
                        exit llpp
                    endif
                endif
             enddo llpp
             dist=0.
             if(pt(1)>0.and.pt(2)>0) &
                 dist=sqrt((vX(pt(1))-vX(pt(2)))**2.0+(vY(pt(2))-vY(pt(1)))**2.0)
             tri_voronoi_edge(cell,n)=dist
             tri_distance(cell,n)=sqrt((tX(cell)-tX(id))**2.0+(tY(cell)-tY(id))**2.0)
          enddo
       endif
    enddo

    return

  end subroutine delaunay_voronoi_duality
  ! =====================================================================================
  subroutine envelope(x,y,n,vertex,nvert)

    integer :: n,vertex(n),nvert,iwk(n)
    integer :: next(20),i,i1,i2,j,jp1,jp2,i2save,i3,i2next

    real(kind=8) :: x(n),y(n),xmax,xmin,ymax,ymin,dist,dmax,dmin,x1,y1
    real(kind=8) :: dx,dy,x2,y2,dx1,dx2,dmax1,dmax2,dy1,dy2,temp,zero

    zero = 0.0_8
    x = x
    y = y
    if(n<2) return
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

    iwk(i1)=-1
    iwk(i2)=-1
    nvert=2
    j=1

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

    if(i>0) goto 60
    iwk(i1)=-1
    iwk(i2)=-1

    goto 40

  end subroutine envelope
  ! =====================================================================================
  function find_intersection_between_segments(pt1,pt2,pt3,pt4) result(solution)

      logical :: solution
      real(kind=8),dimension(2)::pt1,pt2,pt3,pt4
      real(kind=8),dimension(2)::p13,p43,p21
      real(kind=8)::num_a,num_b,denom,mua,mub

      solution=.false.
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
         if(mua>=0.0_8 .and. mua<=1.0_8)then
            solution = .true.
            return
         else
            return
         endif
      endif

   end function find_intersection_between_segments

end module classfv
