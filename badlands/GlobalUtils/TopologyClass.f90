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
!       Filename:  TopologyClass.f90
!
!    Description:  Implements simple functionalities and store global topology parameters
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module topology

  use strings
  use parameters

  implicit none

  ! Shewchuck Triangle's algorithm generated files
  character(len=128)::TINfile,elemfile,vorofile,rstfolder

  ! Total number of nodes and elements of the user provided regular grid
  integer::nbnodes
  integer::bnbnodes
  integer::nbelm

  ! Users requested minimum delaunay triangle angle and maximum area
  real(kind=8)::del_area,del_angle

  ! Boundary type
  integer,dimension(4)::bounds

  ! Outlet type
  integer::outlet

  ! Delaunay cell centroid
  real(kind=8),dimension(:,:),allocatable::dcentroid
  real(kind=8),dimension(:,:),allocatable::Fdata,Fdata2

  ! Shewchuck Triangle's algorithm delaunay elements/edges and voronoi edges
  integer,dimension(:,:),allocatable::delmt,dedg,vedg,delmt2

  ! Shewchuck Triangle's algorithm grid boundaries, voronoi duplicates and output masks
  integer,dimension(:),allocatable::tbound,mask,elemtmask

  ! Number of user defined refinement regions
  integer::refineNb

  ! Users requested refinement regions parameters
  type regions_refine
     real(kind=8),dimension(2)::xcoord
     real(kind=8),dimension(2)::ycoord
     real(kind=8)::area
  end type regions_refine
  type(regions_refine),dimension(:),allocatable::refine_grid

contains
  ! =====================================================================================
  function cmp_centroid(fce) result(centroid)

    integer::fce,k
    integer,dimension(3)::nids

    real(kind=8),dimension(2)::centroid

    nids=delmt(fce,1:3)
    centroid(1:2)=0.0_8
    do k=1,3
       centroid(1)=centroid(1)+tcoordX(nids(k))
       centroid(2)=centroid(2)+tcoordY(nids(k))
    enddo

    do k=1,2
       centroid(k)=centroid(k)/3.0_8
    enddo

  end function cmp_centroid
  ! =====================================================================================
  subroutine ReadStrings(id,str)

    character(len=50)::str
    character(len=1)::delims

    integer,parameter::StrMax=50,Nmax=10
    character(len=StrMax),dimension(Nmax)::args
    integer::nargs,k,id,stat

    delims='   '

    ! split strings according delimiters
    call parsestg(str,delims,args,nargs)

    do k=1,2
       read(args(k+1),*,iostat=stat) vedg(id,k)
       ! If the point is a duplicate associate the edge point to the proper index
       if(vedg(id,k)>0) vedg(id,k)=mask( vedg(id,k))
    enddo

    return

  end subroutine ReadStrings
  ! =====================================================================================
  subroutine addpath(fname)

    integer::pathlen,flen

    character(len=128)::fname,dname,dummy

    ! for files to be read, they'll be in the session path
    dname=' '
    call noblnk(outdir)
    pathlen=len_trim(outdir)
    dname(1:pathlen)=outdir
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath
  ! =====================================================================================
  subroutine addpath1(fname)

    integer::pathlen,flen


    character(len=128)::fname,dname,dummy

    ! for files to be read, they'll be in the session path
    dname=' '
    call noblnk(outdir1)
    pathlen=len_trim(outdir1)
    dname(1:pathlen)=outdir1
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath1
  ! =====================================================================================
  subroutine addpath2(fname)

    integer::pathlen,flen

    character(len=128)::fname,dname,dummy

    ! for files to be read, they'll be in the session path
    dname=' '
    call noblnk(outdir2)
    pathlen=len_trim(outdir2)
    dname(1:pathlen)=outdir2
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath2
  ! =====================================================================================
  subroutine addpath3(fname)

    integer::pathlen,flen

    character(len=128)::fname,dname,dummy

    ! for files to be read, they'll be in the session path
    dname=' '
    call noblnk(outdir3)
    pathlen=len_trim(outdir3)
    dname(1:pathlen)=outdir3
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath3
  ! =====================================================================================
  subroutine noblnk(string)

    integer::i,j,lg

    character(len=128)::string

    lg=len(string)
    do
       if(lg<=0 .or. string(lg:lg)/=' ') exit
       lg=lg-1
    enddo

    if(lg>0)then
       ! find first non-blank character
       i=1
       do
          if(i>lg .or. (string(i:i)/=' ' .and. string/=' ')) exit
          i=i+1
       enddo
       ! determine end of continuous (non-blank) string
       j=i
       do
          if(j>lg)then
             exit
          elseif(string(j:j)==' ')then
             exit
          elseif(string=='  ')then
             exit
          elseif(j==128)then
             exit
          endif
          j=j+1
       enddo
       ! j points to first blank position or past end of string; adjust to last
       ! non-blank position in string
       j=min(j-1,lg)
       string=string(i:j)
       if(j<len(string)) string(j+1:len(string))=' '
    else
       ! there were only blanks in string
       string=' '
    endif

    return

  end subroutine noblnk
  ! ============================================================================
  subroutine append_str2(stg1,stg2)

    integer::l1,l2
    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2

    return

  end subroutine append_str2
  ! =====================================================================================
  subroutine append_str(stg1,stg2)

    integer :: l1,l2
    character(len=128) :: stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_str
  ! =====================================================================================
  subroutine append_zero(stg1,i)

    integer::l1,l2,i

    character(len=128)::stg1,stg2,str

    l2=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l1=len_trim(stg2)
    str=''
    if(l1==1)then
       str(1:3)='000'
       call append_str(str,stg2)
    elseif(l1==2)then
       str(1:2)='00'
       call append_str(str,stg2)
    elseif(l1==3)then
       str(1:1)='0'
       call append_str(str,stg2)
    elseif(l1==4)then
       call append_str(str,stg2)
    endif
    l1=len_trim(str)
    stg1(l2+1:l2+l1)=str
    call noblnk(stg1)

    return

  end subroutine append_zero
  ! =====================================================================================
  subroutine append_nbreal(stg1,n)

    real(kind=8)::n

    integer::l1,l2

    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    write(stg2,'(f30.4)')n
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=trim(stg2)
    call noblnk(stg1)

    return

  end subroutine append_nbreal
  ! =====================================================================================
  subroutine append_nb(stg1,i)

    integer::l1,l2,i

    character(len=128)::stg1,stg2

    l1 = len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_nb
  ! =====================================================================================
  subroutine append_nb2(stg1,i)

    integer::l1,l2,i

    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+2:l1+l2+1)=stg2

    return

  end subroutine append_nb2
  ! =====================================================================================
end module topology
