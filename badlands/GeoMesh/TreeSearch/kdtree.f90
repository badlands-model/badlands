Module kd_tree

   ! K-D tree routines in Fortran 90 by Matt Kennel.
   ! Original program was written in Sather by Steve Omohundro and
   ! Matt Kennel.  Only the Euclidean metric works so far.


   ! Institute for Nonlinear science, UC San Diego, 1997
   ! kennel@lyapunov.ucsd.edu
   ! .. Parameters ..
   ! you choose this.
   Integer, Parameter :: bucket_size = 10
   ! ..
   ! .. Derived Type Declarations ..
   ! Global information about the tree
   ! pointer to the actual
   ! data array
   ! dimensionality and total # of points
   ! permuted index into the
   ! data, so that
   ! indexes[l..u]
   ! of some bucket represent the indexes of the actual
   ! points in that bucket.
   ! root pointer of the tree
   ! an internal tree node
   ! the dimension to cut
   ! where to cut the dimension
   ! indices of points included in this node,
   ! referring back to indices
   ! child pointers
   ! One of these is created for each search.
   ! best squared distance found so far
   ! best index found so far
   ! Query vector
   ! indexes of best found so far
   ! squared distances found so far
   ! how many best distances we are searching
   ! for
   ! how many have been found so far, i.e. il(1:nfound)
   ! are the only valid indexes.
   ! exclude points within
   ! 'correltime'
   ! of 'centeridx'

   Type :: tree_node
      Integer :: dnum
      Real :: val
      Integer :: l, u
      Type (tree_node), Pointer :: left, right
   End Type tree_node

   Type :: tree_master_record
      Real, Dimension (:,:), Pointer :: the_data
      Integer :: dim, n
      Integer, Dimension (:), Pointer :: indexes
      Type (tree_node), Pointer :: root
   End Type tree_master_record

   Type :: tree_search_record
      Private
      Real :: bsd
      Real, Dimension (:), Pointer :: qv
      Integer, Dimension (:), Pointer :: il
      Real, Dimension (:), Pointer :: dsl
      Integer :: nbst, nfound
      Integer :: centeridx, correltime
      Logical :: linfinity_metric
   End Type tree_search_record
  ! ..
  ! set to true if we're doing linfinity metric
Contains

   Subroutine destroy_tree(tp)
      ! Deallocates all memory for the tree, except input data matrix
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      Call destroy_node(tp%root)
      Deallocate (tp%indexes)
      tp%indexes => NULL() !#Nullify (tp%indexes)
      Return
   end Subroutine destroy_tree
   Recursive Subroutine destroy_node(np)
      ! .. Structure Arguments ..
      Type (tree_node), Pointer :: np
      ! ..
      ! ..
      If (ASSOCIATED(np%left)) Then
         Call destroy_node(np%left)
         Deallocate (np%left)
         np%left => NULL()
      End If
      If (ASSOCIATED(np%right)) Then
         Call destroy_node(np%right)
         Deallocate (np%right)
         np%right => NULL()
      End If
      Return
   End Subroutine destroy_node


   Function create_tree(input_data) Result (master_record)
      ! create the actual tree structure, given an input array of data.
      ! Arguments
      ! .. Function Return Value ..
      Type (tree_master_record), Pointer :: master_record
      ! ..
      ! .. Array Arguments ..
      Real, Target :: input_data(:,:)
      ! ..
      ! ..
      Allocate (master_record)
      master_record%the_data => input_data
      ! pointer assignment
      master_record%n = SIZE(input_data,1)
      master_record%dim = SIZE(input_data,2)
      Print *, 'Creating KD tree with N = ', master_record%n, &
      ' and dim = ', master_record%dim
      Call build_tree(master_record)

   Contains

      Subroutine build_tree(tp)
         ! .. Structure Arguments ..
         Type (tree_master_record), Pointer :: tp
         ! ..
         ! .. Local Scalars ..
         Integer :: j
         ! ..
         Allocate (tp%indexes(tp%n))
         Do j = 1, tp%n
            tp%indexes(j) = j
         End Do
         tp%root => build_tree_for_range(tp,1,tp%n)
      End Subroutine build_tree

      Recursive Function build_tree_for_range(tp,l,u) Result (res)
         ! .. Function Return Value ..
         Type (tree_node), Pointer :: res
         ! ..
         ! .. Structure Arguments ..
         Type (tree_master_record), Pointer :: tp
         ! ..
         ! .. Scalar Arguments ..
         Integer, Intent (In) :: l, u
         ! ..
         ! .. Local Scalars ..
         Integer :: c, m
         ! ..
         If ((u-l)<=bucket_size) Then
            Allocate (res)
            res%dnum = 0
            res%val = 0.0
            res%l = l
            res%u = u
            Nullify (res%left,res%right)
         Else
            Allocate (res)
            c = most_spread_coordinate(tp,l,u)
            m = (l+u)/2
            Call select_on_coordinate(tp,c,m,l,u)
            ! moves indexes around
            res%dnum = c
            res%val = tp%the_data(tp%indexes(m),c)
            res%l = l
            res%u = u
            res%left => build_tree_for_range(tp,l,m)
            res%right => build_tree_for_range(tp,m+1,u)
         End If
      End Function build_tree_for_range

      Subroutine select_on_coordinate(tp,c,k,li,ui)
         ! Move elts of ind around between l and u, so that the kth
         ! element
         ! is >= those below, <= those above, in the coordinate c.
         ! .. Structure Arguments ..
         Type (tree_master_record), Pointer :: tp
         ! ..
         ! .. Scalar Arguments ..
         Integer, Intent (In) :: c, k, li, ui
         ! ..
         ! .. Local Scalars ..
         Integer :: i, l, m, s, t, u
         ! ..
         ! .. Local Arrays ..
         Real, Pointer :: v(:,:)
         Integer, Pointer :: ind(:)
         ! ..
         v => tp%the_data
         ind => tp%indexes
         l = li
         u = ui
         Do While (l<u)
            t = ind(l)
            m = l
            Do i = l + 1, u
               If (v(ind(i),c)<v(t,c)) Then
                  m = m + 1
                  s = ind(m)
                  ind(m) = ind(i)
                  ind(i) = s
               End If
            End Do
            s = ind(l)
            ind(l) = ind(m)
            ind(m) = s
            If (m<=k) l = m + 1
            If (m>=k) u = m - 1
         End Do
      End Subroutine select_on_coordinate

      Function most_spread_coordinate(tp,l,u) Result (res)
         ! Of indices in l..u find the axis which has the largest spread,
         ! and
         ! return its index
         ! .. Function Return Value ..
         Integer :: res
         ! ..
         ! .. Structure Arguments ..
         Type (tree_master_record), Pointer :: tp
         ! ..
         ! .. Scalar Arguments ..
         Integer, Intent (In) :: l, u
         ! ..
         ! .. Local Scalars ..
         Real :: bsp, sp
         Integer :: i = 0
         ! ..
         res = 0
         bsp = -1.0
         Do i = 1, tp%dim
            sp = spread_in_coordinate(tp,i,l,u)
            If (sp>bsp) Then
               res = i
               bsp = sp
            End If
         End Do
      End Function most_spread_coordinate

      Function spread_in_coordinate(tp,c,l,u) Result (res)
         ! the spread in coordinate 'c', between l and u
         ! for easier local access
         ! ibid
         ! .. Function Return Value ..
         Real :: res
         ! ..
         ! .. Structure Arguments ..
         Type (tree_master_record), Pointer :: tp
         ! ..
         ! .. Scalar Arguments ..
         Integer, Intent (In) :: c, l, u
         ! ..
         ! .. Local Scalars ..
         Real :: last, lmax, lmin, max, min, t
         Integer :: i
         ! ..
         ! .. Local Arrays ..
         Real, Pointer :: v(:,:)
         Integer, Pointer :: ind(:)
         ! ..
         v => tp%the_data
         ind => tp%indexes
         min = v(ind(l),c)
         max = min
         Do i = l + 2, u, 2
            lmin = v(ind(i-1),c)
            lmax = v(ind(i),c)
            If (lmin>lmax) Then
               t = lmin
               lmin = lmax
               lmax = t
            End If
            If (min>lmin) min = lmin
            If (max<lmax) max = lmax
         End Do
         If (i==u+1) Then
            last = v(ind(u),c)
            If (min>last) min = last
            If (max<last) max = last
         End If
         res = max - min
        ! write *, "Returning spread in coordinate with res = ", res
      End Function spread_in_coordinate
   End Function create_tree

   ! Search routines:

   ! * n_nearest_to(tp, qv, n, indexes, distances)

   ! Find the 'n' vectors in the tree nearest to 'qv' in euclidean norm
   ! returning their indexes and distances in 'indexes' and 'distances'
   ! arrays already allocated passed to the subroutine.
   Subroutine n_nearest_to(tp,qv,n,indexes,distances)
      ! find the 'n' nearest neighbors to 'qv', returning
      ! their indexes in indexes and squared Euclidean distances in
      ! 'distances'
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: n
      ! ..
      ! .. Array Arguments ..
      Real, Target :: distances(n)
      Real, Target, Intent (In) :: qv(:)
      Integer, Target :: indexes(n)
      ! ..
      ! .. Local Structures ..
      Type (tree_search_record), Pointer :: psr
      Type (tree_search_record), Target :: sr
      ! ..
      ! ..
      ! the largest real number
      sr%bsd = HUGE(1.0)
      sr%qv => qv
      sr%nbst = n
      sr%nfound = 0
      sr%centeridx = -1
      sr%correltime = 0
      sr%linfinity_metric = .False. ! change if you want.
      sr%dsl => distances
      sr%il => indexes
      sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
      sr%il = -1               ! set to invalid indexes
      psr => sr                ! in C this would be psr = &sr
      Call n_search(tp,psr,tp%root)
      Return
   End Subroutine n_nearest_to

   ! Another main routine you call from your user program


   Subroutine n_nearest_to_around_point(tp,idxin,correltime,n,indexes, &
   distances)
      ! find the 'n' nearest neighbors to point 'idxin', returning
      ! their indexes in indexes and squared Euclidean distances in
      ! 'distances'
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer, Intent (In) :: correltime, idxin, n
      ! ..
      ! .. Array Arguments ..
      Real, Target :: distances(n)
      Integer, Target :: indexes(n)
      ! ..
      ! .. Local Structures ..
      Type (tree_search_record), Pointer :: psr
      Type (tree_search_record), Target :: sr
      ! ..
      ! .. Local Arrays ..
      Real, Allocatable, Target :: qv(:)
      ! ..
      ! ..
      Allocate (qv(tp%dim))
      qv = tp%the_data(idxin,:) ! copy the vector
      sr%bsd = HUGE(1.0)       ! the largest real number
      sr%qv => qv
      sr%nbst = n
      sr%nfound = 0
      sr%centeridx = idxin
      sr%correltime = correltime
      sr%linfinity_metric = .False.
      sr%dsl => distances
      sr%il => indexes
      sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
      sr%il = -1               ! set to invalid indexes
      psr => sr                ! in C this would be psr = &sr
      Call n_search(tp,psr,tp%root)
      Deallocate (qv)
      Return
   End Subroutine n_nearest_to_around_point

   Function bounded_square_distance(tp,qv,ii,bound) Result (res)
      ! distance between v[i,*] and qv[*] if less than or equal to bound,
      ! -1.0 if greater than than bound.
      ! .. Function Return Value ..
      Real :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Real :: bound
      Integer :: ii
      ! ..
      ! .. Array Arguments ..
      Real :: qv(:)
      ! ..
      ! .. Local Scalars ..
      Real :: tmp
      Integer :: i
      ! ..
      res = 0.0
      Do i = 1, tp%dim
         tmp = tp%the_data(ii,i) - qv(i)
         res = res + tmp*tmp
         If (res>bound) Then
            res = -1.0
            Return
         End If
      End Do
   End Function bounded_square_distance

   Recursive Subroutine n_search(tp,sr,node)
      ! This is the innermost core routine of the kd-tree search.
      ! it is thus not programmed in quite the clearest style, but is
      ! designed for speed.  -mbk
      ! .. Structure Arguments ..
      Type (tree_node), Pointer :: node
      Type (tree_search_record), Pointer :: sr
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Local Scalars ..
      Real :: dp, sd, sdp, tmp
      Integer :: centeridx, i, ii, j, jmax, k
      Logical :: condition, not_fully_sized
      ! ..
      ! .. Local Arrays ..
      Real, Pointer :: qv(:)
      Integer, Pointer :: ind(:)
      ! ..
      ! ..
      If ( .Not. (ASSOCIATED(node%left)) .And. ( .Not. ASSOCIATED( &
      node%right))) Then
         ! we are on a terminal node
         ind => tp%indexes     ! save for easy access
         qv => sr%qv
         centeridx = sr%centeridx
         ! search through terminal bucket.
         Do i = node%l, node%u
            ii = ind(i)
            condition = .False.
            If (centeridx<0) Then
               condition = .True.
            Else
               condition = (ABS(ii-sr%centeridx)>=sr%correltime)
            End If
            If (condition) Then
               ! sd = bounded_square_distance(tp,qv,ii,sr%bsd)
               ! replace with call to bounded_square distance with inline
               ! code, an
               ! a goto.  Tested to be significantly faster.   SPECIFIC
               ! FOR
               ! the EUCLIDEAN METRIC ONLY! BEWARE!

               sd = 0.0
               Do k = 1, tp%dim
                  tmp = tp%the_data(ii,k) - qv(k)
                  sd = sd + tmp*tmp
                  If (sd>sr%bsd) Then
                     Go To 100
                  End If
               End Do
               ! we only consider it if it is better than the 'best' on
               ! the list so far.
               ! if we get here
               ! we know sd is < bsd, and bsd is on the end of the list
               not_fully_sized = (sr%nfound<sr%nbst)
               If (not_fully_sized) Then
                  jmax = sr%nfound
                  sr%nfound = sr%nfound + 1
               Else
                  jmax = sr%nbst - 1
               End If
               ! add it to the list
               Do j = jmax, 1, -1
                  If (sd>=sr%dsl(j)) Exit ! we hit insertion location
                  sr%il(j+1) = sr%il(j)
                  sr%dsl(j+1) = sr%dsl(j)
               End Do
               sr%il(j+1) = ii
               sr%dsl(j+1) = sd
               If ( .Not. not_fully_sized) Then
                  sr%bsd = sr%dsl(sr%nbst)
               End If
100         Continue        ! we jumped  here if we had a bad point.
         End If
      End Do
   Else
      ! we are not on a terminal node

      ! Alrighty, this section is essentially the content of the
      ! the Sproul method for searching a Kd-tree, in other words
      ! the second nested "if" statements in the two halves below
      ! and when they get activated.
      dp = sr%qv(node%dnum) - node%val
      If ( .Not. sr%linfinity_metric) Then
         sdp = dp*dp        ! Euclidean
      Else
         sdp = ABS(dp)      ! Linfinity
      End If
      If (dp<0.0) Then
         Call n_search(tp,sr,node%left)
         If (sdp<sr%bsd) Call n_search(tp,sr,node%right)
         ! if the distance projected to the 'wall boundary' is less
         ! than the radius of the ball we must consider, then perform
         ! search on the 'other side' as well.
      Else
         Call n_search(tp,sr,node%right)
         If (sdp<sr%bsd) Call n_search(tp,sr,node%left)
      End If
   End If
End Subroutine n_search

Subroutine n_nearest_to_brute_force(tp,qv,n,indexes,distances)
   ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
   ! only use this subroutine for testing, as it is SLOW!  The
   ! whole point of a k-d tree is to avoid doing what this subroutine
   ! does.
   ! .. Structure Arguments ..
   Type (tree_master_record), Pointer :: tp
   ! ..
   ! .. Scalar Arguments ..
   Integer, Intent (In) :: n
   ! ..
   ! .. Array Arguments ..
   Real, Target :: distances(n)
   Real, Intent (In) :: qv(:)
   Integer, Target :: indexes(n)
   ! ..
   ! .. Local Scalars ..
   Integer :: i, j, k
   ! ..
   ! .. Local Arrays ..
   Real, Allocatable :: all_distances(:)
   ! ..
   ! ..
   Allocate (all_distances(tp%n))
   Do i = 1, tp%n
      all_distances(i) = bounded_square_distance(tp,qv,i,HUGE(1.0))
   End Do
   ! now find 'n' smallest distances
   Do i = 1, n
      distances(i) = HUGE(1.0)
      indexes(i) = -1
   End Do
   Do i = 1, tp%n
      If (all_distances(i)<distances(n)) Then
         ! insert it somewhere on the list
         Do j = 1, n
            If (all_distances(i)<distances(j)) Exit
         End Do
         ! now we know 'j'
         Do k = n - 1, j, -1
            distances(k+1) = distances(k)
            indexes(k+1) = indexes(k)
         End Do
         distances(j) = all_distances(i)
         indexes(j) = i
      End If
   End Do
   Deallocate (all_distances)
End Subroutine n_nearest_to_brute_force
End Module kd_tree
