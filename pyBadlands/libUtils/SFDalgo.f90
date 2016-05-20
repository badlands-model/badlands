!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! This module computes the Single Flow Direction for any given surface.

module sfdcompute

  implicit none

contains

    subroutine directions(pyElev,pyZ,pyNgbs,pyEdge,pyDist,pyGIDs,sealimit,pyBase, &
        pyRcv,pyMaxh,pyMaxDep,pyDiff,pylocalNb,pyglobalNb)

        integer :: pylocalNb
        integer :: pyglobalNb
        real(kind=8),intent(in) :: sealimit
        integer,dimension(pylocalNb),intent(in) :: pyGIDs
        integer,dimension(pyglobalNb,20),intent(in) :: pyNgbs
        real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
        real(kind=8),dimension(pyglobalNb),intent(in) :: pyElev
        real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
        real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist

        integer,intent(out) :: pyBase(pyglobalNb)
        integer,intent(out) :: pyRcv(pyglobalNb)
        real(kind=8),intent(out) :: pyMaxh(pyglobalNb)
        real(kind=8),intent(out) :: pyDiff(pyglobalNb)
        real(kind=8),intent(out) :: pyMaxDep(pyglobalNb)

        integer :: p,k,lowestID,gid
        real(kind=8) :: diffH,diffD,dh

        pyBase = -1
        pyRcv = -1
        pyMaxh = 1.e6
        pyMaxDep = 0.
        pyDiff = 0.
        do k = 1, pylocalNb
            gid = pyGIDs(k)+1
            lowestID = gid
            diffH = 1.e6
            diffD = 0.
            p = 1
            do while(pyNgbs(gid,p) >=0 )
                if(pyElev(pyNgbs(gid,p)+1) < pyElev(lowestID))then
                    lowestID = pyNgbs(gid,p)+1
                endif
                dh = pyZ(pyNgbs(gid,p)+1)-pyZ(gid)
                if(dh >= 0.) diffH = min(dh, diffH)
                diffD = max(dh, diffD)
                pyDiff(gid) = pyDiff(gid) + pyEdge(gid,p)*dh/pyDist(gid,p)
                p = p+1
            enddo
            pyRcv(gid) = lowestID-1
            if( pyZ(gid) < sealimit ) pyRcv(gid) = gid-1
            if( gid == pyRcv(gid)+1 .and. pyZ(gid)+diffD > sealimit) pyBase(gid) = gid-1
            if( diffH > 9.99e5 ) diffH = 0.
            pyMaxh(gid) = diffH
            pyMaxDep(gid) = diffD
        enddo

        return

    end subroutine directions

    subroutine directions_base(pyZ,pyNgbs,pyEdge,pyDist,pyGIDs,sealimit,pyBase, &
        pyRcv,pyDiff,pylocalNb,pyglobalNb)

        integer :: pylocalNb
        integer :: pyglobalNb
        real(kind=8),intent(in) :: sealimit
        integer,dimension(pylocalNb),intent(in) :: pyGIDs
        integer,dimension(pyglobalNb,20),intent(in) :: pyNgbs
        real(kind=8),dimension(pyglobalNb),intent(in) :: pyZ
        real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyEdge
        real(kind=8),dimension(pyglobalNb,20),intent(in) :: pyDist

        integer,intent(out) :: pyBase(pyglobalNb)
        integer,intent(out) :: pyRcv(pyglobalNb)
        real(kind=8),intent(out) :: pyDiff(pyglobalNb)

        integer :: p,k,lowestID,gid
        real(kind=8) :: dh, diffD

        pyBase = -1
        pyRcv = -1
        pyDiff = 0.
        do k = 1, pylocalNb
            gid = pyGIDs(k)+1
            lowestID = gid
            diffD = 0.
            p = 1
            do while(pyNgbs(gid,p) >=0 )
                if(pyZ(pyNgbs(gid,p)+1) < pyZ(lowestID))then
                    lowestID = pyNgbs(gid,p)+1
                endif
                dh = pyZ(pyNgbs(gid,p)+1)-pyZ(gid)
                diffD = max(dh, diffD)
                pyDiff(gid) = pyDiff(gid) + pyEdge(gid,p)*dh/pyDist(gid,p)
                p = p+1
            enddo
            pyRcv(gid) = lowestID-1
            if( pyZ(gid) < sealimit ) pyRcv(gid) = gid-1
            if( gid == pyRcv(gid)+1 .and. pyZ(gid)+diffD > sealimit) pyBase(gid) = gid-1
        enddo

        return

    end subroutine directions_base

end module sfdcompute
