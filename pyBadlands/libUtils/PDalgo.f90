!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module implements Planchon & Darboux depression filling algorithm

module pdcompute

  implicit none

contains

  subroutine filling(elevation,pyNgbs,fillTH,epsilon,pybounds,sealimit,demH,pydnodes)

      integer :: pydnodes
      integer,intent(in) :: pybounds
      real(kind=8),intent(in) :: sealimit
      real(kind=8),intent(in) :: fillTH
      real(kind=8),intent(in) :: epsilon
      real(kind=8),intent(in) :: elevation(pydnodes)
      integer,intent(in) :: pyNgbs(pydnodes,20)
      real(kind=8),intent(out) :: demH(pydnodes)

      logical :: flag
      integer :: p, k

      flag = .true.
      demH = 1.e6
      demH(1:pybounds) = elevation(1:pybounds)
      do k=pybounds+1,pydnodes
          if( demH(k) > sealimit)then
              demH(k) = 1.e6
          else
              demH(k) = elevation(k)
          endif
      enddo

      do while(flag)
        flag=.false.
        do k=1,pydnodes
          if( demH(k) > elevation(k) )then
            loop: do p = 1, 20
                if( pyNgbs(k,p) < 0) exit loop
                if( elevation(k) >= demH(pyNgbs(k,p)+1) + epsilon )then
                    demH(k) = elevation(k)
                else
                    if( demH(k) > demH(pyNgbs(k,p)+1) + epsilon )then
                        demH(k) = demH(pyNgbs(k,p)+1) + epsilon
                        if(demH(k) - elevation(k) > fillTH)then
                            demH(k) = elevation(k) + fillTH
                        else
                            flag=.true.
                        endif
                    endif
                endif
            enddo loop
          endif
        enddo
      enddo

      return

  end subroutine filling

  ! THIS IS NOT USED
  subroutine basin_filling(pyGIDS,pyLIDs,elevation,pyNgbs,epsilon,pybounds,demH,pylNodesNb,pygNodesNb)

      integer :: pylNodesNb
      integer :: pygNodesNb
      integer,intent(in) :: pybounds
      real(kind=8),intent(in) :: epsilon
      real(kind=8),intent(in) :: elevation(pygNodesNb)
      integer,intent(in) :: pyGIDS(pylNodesNb)
      integer,intent(in) :: pyLIDS(pygNodesNb)
      integer,intent(in) :: pyNgbs(pylNodesNb,20)
      real(kind=8),intent(out) :: demH(pylNodesNb)

      logical :: flag
      integer :: p, k, gid, lid

      flag = .true.
      demH = 1.e6
      do k = 1, pybounds
        gid = pyGIDS(k) + 1
        demH(k) = elevation(gid)
      enddo

      do while(flag)
        flag=.false.
        do k=1,pylNodesNb
          gid = pyGIDS(k) + 1
          if( demH(k) > elevation(gid) )then
            loop: do p = 1, 20
              if( pyNgbs(k,p) < 0) exit loop
              lid = pyLIDS(pyNgbs(k,p)+1) + 1
              if( elevation(gid) >= demH(lid) + epsilon )then
                demH(k) = elevation(gid)
              else
                if( demH(k) > demH(lid) + epsilon )then
                  demH(k) = demH(lid) + epsilon
                  flag=.true.
                endif
              endif
            enddo loop
          endif
        enddo
      enddo

      return

  end subroutine basin_filling

  ! THIS IS NOT USED
  subroutine fill_recursive(depFill, pyArea, pyElev, distH, pybaseNb)

    integer :: pybaseNb
    real(kind=8),intent(in) :: depFill
    real(kind=8),intent(in) :: pyArea(pybaseNb)
    real(kind=8),intent(in) :: pyElev(pybaseNb)

    real(kind=8),intent(out) :: distH(pybaseNb)

    integer :: id, stpID
    real(kind=8) :: dh, cumFill, area, fillh, distVol, flatVol

    id = 1
    area = 0.
    cumFill = 0.
    distVol = 0.01 * depFill
    flatVol = depFill !- distVol

    fillh = pyElev(1)
    lp: do id = 1, pybaseNb-1

      stpID = id

      dh = pyElev(id+1) - pyElev(id)
      area = area + pyArea(id)

      if( cumFill + area * dh < flatVol )then
        fillh = fillh + dh
        cumFill = cumFill + area * dh
      else
        fillh = fillh + ( flatVol - cumFill )/area
        exit lp
      endif

    enddo lp

    do id = 1, pybaseNb
      if(fillh - pyElev(id) > 0)then
        distH(id) = fillh - pyElev(id)
      else
        distH(id) = 0.
      endif
    enddo

    return

    dh = 0.01
    distH = fillh
    lp2: do while(distVol > 0.)
      do id = stpID, 2
        distH(id) = distH(id) + dh * (id/stpID)
        distVol = distVol - dh * pyArea(id) * (id/stpID)
        if(distVol < 0.) exit lp2
      enddo
    enddo lp2

    return

  end subroutine fill_recursive

end module pdcompute
