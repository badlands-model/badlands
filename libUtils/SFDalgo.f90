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
    
    subroutine directions(pyElev,pyNgbs,pyGIDs,pyBase,pyRcv,pylocalNb,pyglobalNb)

        integer :: pylocalNb
        integer :: pyglobalNb
        
        integer,dimension(pylocalNb),intent(in) :: pyGIDs
        integer,dimension(pylocalNb,20),intent(in) :: pyNgbs
        real,dimension(pyglobalNb),intent(in) :: pyElev
        
        integer,intent(out) :: pyBase(pyglobalNb)
        integer,intent(out) :: pyRcv(pyglobalNb)
        
        integer :: p,k,lowestID,gid
        real :: minElev
        
        pyBase = -1
        pyRcv = -1
        do k = 1, pylocalNb
            gid = pyGIDs(k)+1
            lowestID = gid
            minElev = 1.e7
            p = 1
            do while(pyNgbs(k,p)>=0)
                if(pyElev(pyNgbs(k,p)+1) < pyElev(lowestID))then 
                    lowestID = pyNgbs(k,p)+1
                endif
                p = p+1
            enddo
            pyRcv(gid) = lowestID-1
            if( gid == pyRcv(gid)+1 ) pyBase(gid) = gid-1
        enddo
        
        return

    end subroutine directions

end module sfdcompute