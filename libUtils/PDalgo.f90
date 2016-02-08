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

  use pdstack
  
  implicit none

contains

  subroutine filling(elevation,pyNgbs,fillTH,epsilon,pybounds,demH,pydnodes)

      integer :: pydnodes 
      integer,intent(in) :: pybounds
      real(kind=8),intent(in) :: fillTH
      real(kind=8),intent(in) :: epsilon
      real(kind=8),intent(in) :: elevation(pydnodes)
      integer,intent(in) :: pyNgbs(pydnodes,20)
      real(kind=8),intent(out) :: demH(pydnodes)
       
      logical::flag
      integer::p,k,l
      real(kind=8)::minH
      
      stack1%size = 0
      stack2%size = 0
      demH(1:pybounds) = elevation(1:pybounds)
      do k=pybounds+1,pydnodes
          demH(k) = 1.e6
          call push(stack1,k)
      enddo
      
      flag = .true.
      do while( flag )
          flag = .false.
          do l = 1,stack1%size
              k = stack1%data(l)
              if( demH(k) > elevation(k) )then
                  minH = demH(pyNgbs(k,1)+1)
                  p = 2
                  do while(pyNgbs(k,p) >= 0) 
                      minH = min(demH(pyNgbs(k,p)+1),minH)
                      p = p+1
                  enddo
                  if (elevation(k) >= minH + epsilon )then
                       demH(k) = elevation(k)
                  else
                      if( demH(k) > minH + epsilon )then
                          demH(k) = minH + epsilon
                          if( demH(k) - elevation(k) >= fillTH )then
                              demH(k) = elevation(k) + fillTH
                          else
                              call push(stack2,k)
                              flag = .true.
                          endif
                      else
                          call push(stack2,k)
                      endif
                  endif
              endif
          enddo
          if(flag)then 
              stack1%size = stack2%size
              stack1%data(1:stack1%size ) = stack2%data(1:stack2%size)
              stack2%size = 0
          endif
       enddo
       
       return
       
  end subroutine filling

end module pdcompute