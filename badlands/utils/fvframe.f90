!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!



subroutine definetin( coords, cells_nodes, cells_edges, edges_nodes, &
                      circumcenter, ngbid, vor_edges, edge_length, maxngbhs, n, nb, m)
!*****************************************************************************
! Compute for a specific triangulation the characteristics of each node and
! associated voronoi for finite volume discretizations

  implicit none

  integer :: m, n, nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: cells_edges(n,3)
  integer, intent(in) :: edges_nodes(m, 2)

  real( kind=8 ), intent(in) :: coords(nb,3)
  real( kind=8 ), intent(in) :: circumcenter(3,n)


  integer, intent(out) :: ngbid(nb, 20)
  real( kind=8 ), intent(out) :: edge_length(nb, 20)
  real( kind=8 ), intent(out) :: vor_edges(nb, 20)
  integer, intent(out) :: maxNgbhs

  integer :: fvnnb(nb)
  integer :: i, n1, n2, k, l, p, eid, cid, e, id
  integer :: nid(2), nc(3), edge(nb, 8)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 20)

  real( kind=8 ) :: coords0(3), coordsID(3)
  real( kind=8 ) :: midpoint(3), dist

  logical :: inside

  cell_ids = -1
  edge = -1
  fvnnb = 0
  ngbid = -1
  edge_length = 0.
  vor_edges = 0.
  maxNgbhs = 0

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells_nodes(i,1:3)+1
    do p = 1, 3
      inside = .False.
      lp: do k = 1, 8
        if( cell_ids(nc(p),k) == i-1 )then
          exit lp
        elseif( cell_ids(nc(p),k) == -1 )then
          inside = .True.
          exit lp
        endif
      enddo lp
      if( inside )then
        cell_ids(nc(p),k)  = i-1
      endif
    enddo
  enddo

  ! Find all edges connected to a given vertice
  do i = 1, m
    n1 = edges_nodes(i,1)+1
    n2 = edges_nodes(i,2)+1
    inside = .False.
    lp0: do k = 1, 8
      if(edge(n1,k) == i-1)then
        exit lp0
      elseif(edge(n1,k) == -1)then
        inside = .True.
        exit lp0
      endif
    enddo lp0
    if( inside )then
      edge(n1,k)  = i-1
      fvnnb(n1) = fvnnb(n1) + 1
      maxNgbhs = max(maxNgbhs,fvnnb(n1))
    endif
    inside = .False.
    lp1: do k = 1, 8
      if(edge(n2,k) == i-1)then
        exit lp1
      elseif(edge(n2,k) == -1)then
        inside = .True.
        exit lp1
      endif
    enddo lp1
    if( inside )then
      edge(n2,k)  = i-1
      fvnnb(n2) = fvnnb(n2) + 1
      maxNgbhs = max(maxNgbhs,fvnnb(n2))
    endif
  enddo

  do k = 1, nb
    ! Get triangulation edge lengths
    coords0 = coords(k,1:3)
    l = 0
    do eid = 1, fvnnb(k)
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        ngbid(k,l) = nid(2)
        coordsID = coords(nid(2)+1,1:3)
      else
        l = l + 1
        ngbid(k,l) = nid(1)
        coordsID = coords(nid(1)+1,1:3)
      endif
      call euclid( coords0, coordsID, edge_length(k,l) )
    enddo

    ! Get voronoi edge lengths
    lp2: do cid = 1, 20
      if( cell_ids(k,cid) == -1 ) exit lp2
      edgeNb(1:3) = cells_edges( cell_ids(k,cid)+1,1:3 )
      do e = 1, 3
        edges(e,1:2) = edges_nodes(edgeNb(e)+1,1:2)
        if( k-1 == edges(e,1) .or. k-1 == edges(e,2))then
          midpoint(1:3) = 0.5 * (coords(edges(e,1)+1,1:3)+coords(edges(e,2)+1,1:3))
          id = -1
          if( edges(e,1) == k-1 )then
            lp3: do i = 1, fvnnb(k)
              if(ngbid(k,i) == edges(e,2))then
                id = i
                exit lp3
              endif
            enddo lp3
          else
            lp4: do i = 1, fvnnb(k)
              if(ngbid(k,i) == edges(e,1))then
                id = i
                exit lp4
              endif
            enddo lp4
          endif
          if( id > -1)then
            call euclid( midpoint(1:3), circumcenter(1:3,cell_ids(k,cid)+1),  dist)
            vor_edges(k,id) = vor_edges(k,id) + dist
          endif
        endif
      enddo
    enddo lp2

  enddo

end subroutine definetin

subroutine euclid( p1, p2, norm)
!*****************************************************************************
! Computes the Euclidean vector norm between 2 points
! along the 3 dimensions

  implicit none

  real( kind=8 ), intent(in) :: p1(3)
  real( kind=8 ), intent(in) :: p2(3)
  real( kind=8 ), intent(out) :: norm
  real( kind=8 ) :: vec(3)

  vec = p2 - p1
  norm = norm2(vec)

  return

end subroutine euclid

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
  if(1>2) return
end subroutine build
