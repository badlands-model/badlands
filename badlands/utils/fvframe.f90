!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Finite Volume discretisation implementation for unstructured grid.
subroutine build(pyOids,pyGids,pytX,pytY,pytEdge,pytElmt,pyvX,pyvY,pyvEdge,pyVarea,pyNgbs, &
  pyVlenght,pyDlenght,pymaxNgbh,pygnodes,pyonodes,pydnodes,pydedges,pydelems,pyvnodes,pyvedges)

  use classfv
  implicit none

  integer :: pyonodes
  integer :: pygnodes
  integer :: pydnodes
  integer :: pyvnodes
  integer :: pydelems
  integer :: pydedges
  integer :: pyvedges
  integer,intent(in) :: pyOids(pyonodes)
  integer,intent(in) :: pyGids(pygnodes)
  !integer,intent(in) :: pyParts(pydnodes)
  real(kind=8),intent(in) :: pytX(pydnodes)
  real(kind=8),intent(in) :: pytY(pydnodes)
  integer,intent(in) :: pytElmt(pydelems,3)
  integer,intent(in) :: pytEdge(pydedges,2)
  real(kind=8),intent(in) :: pyvX(pyvnodes)
  real(kind=8),intent(in) :: pyvY(pyvnodes)
  integer,intent(in) :: pyvEdge(pyvedges,2)
  ! To compute within the module
  real(kind=8),intent(out) :: pyVarea(pydnodes)
  integer,intent(out) :: pyNgbs(pydnodes,20)
  real(kind=8),intent(out) :: pyVlenght(pydnodes,20)
  real(kind=8),intent(out) :: pyDlenght(pydnodes,20)
  integer,intent(out) :: pymaxNgbh

  onodes = pyonodes
  gnodes = pygnodes
  dnodes = pydnodes
  vnodes = pyvnodes
  delems = pydelems
  dedges = pydedges
  vedges = pyvedges
  pymaxNgbh = 0

  if(allocated(gIDs)) deallocate(gIDs)
  allocate(gIDs(gnodes))
  gIDs = pyGids

  if(allocated(oIDs)) deallocate(oIDs)
  allocate(oIDs(onodes))
  oIDs = pyOids

  if(allocated(tX)) deallocate(tX,tY)
  if(allocated(vX)) deallocate(vX,vY)
  allocate(tX(vnodes),tY(vnodes))
  allocate(vX(vnodes),vY(vnodes))
  tX = pytX
  tY = pytY
  vX = pyvX
  vY = pyvY

  if(allocated(tElmt)) deallocate(tElmt)
  if(allocated(tEdge)) deallocate(tEdge)
  if(allocated(vEdge)) deallocate(vEdge)
  allocate(tElmt(delems,3),tEdge(dedges,2),vEdge(vedges,2))
  tElmt = pytElmt
  tEdge = pytEdge
  vEdge = pyvEdge

  call get_data
  call delaunay_voronoi_duality(.true.)
  call build_export_arrays_o(pyVarea,pyNgbs,pyVlenght,pyDlenght,pymaxNgbh)

end subroutine build
