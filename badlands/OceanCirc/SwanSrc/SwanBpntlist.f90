subroutine SwanBpntlist
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2011  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.80: Marcel Zijlema
!   40.92: Marcel Zijlema
!   41.14: Nico Booij
!   41.20: Casey Dietrich
!
!   Updates
!
!   40.80, April 2008: New subroutine
!   40.92,  June 2008: changes with respect to boundary polygons
!   41.14,  July 2010: boundary segments added as output curve
!   41.20, March 2010: extension to tightly coupled ADCIRC+SWAN model
!
!   Purpose
!
!   Makes list of boundary vertices in ascending order
!   - counterclockwise in case of sea/mainland boundaries
!   - clockwise in case of island boundaries
!
!   Method
!
!   The grid contains a number of boundary polygons
!   They are by definition closed
!   The first boundary polygon refers to sea/mainland boundary and the other polygons refers to island boundaries
!
!   The vertices which define the sea/mainland boundary are inserted in the counterclockwise direction
!   The vertices which define the island boundary are inserted in the clockwise direction
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
    use OUTP_DATA                  ! 41.14
!
    implicit none
!
!   Local variables
!
    integer                            :: icell      ! cell index
    integer, save                      :: ient = 0   ! number of entries in this subroutine
    integer                            :: iface      ! face index
    integer                            :: istat      ! indicate status of allocation
    integer                            :: j          ! loop counter
    integer                            :: k          ! counter
    integer, dimension(1)              :: kx         ! location of minimum value in array of x-coordinates of boundary vertices
    integer, dimension(1)              :: ky         ! location of minimum value in array of y-coordinates of boundary vertices
    integer                            :: m          ! loop counter
    integer                            :: maxnbp     ! maximum number of boundary vertices in set of polygons
    integer                            :: nbptot     ! total number of boundary vertices
    integer                            :: nptemp     ! auxiliary integer to store number of points temporarily
    integer, dimension(3)              :: v          ! vertices in present cell
    integer                            :: v1         ! first vertex of present face
    integer                            :: v2         ! second vertex of present face
    integer                            :: vc         ! considered vertex
    integer                            :: vcf        ! first considered vertex of a boundary polygon
    integer                            :: vn         ! next vertex with respect to considered vertex (counterclockwise)
    integer                            :: MIP        ! number of points on a output curve               41.14
    integer                            :: vmk        ! marker value of a boundary point                 41.14
    integer                            :: JJ         ! counter of points on a curve                     41.14
    integer                            :: VM         ! index of a boundary part                         41.14
    integer                            :: VMMAX      ! highest value of VM                              41.14
    integer                            :: JBG        ! index of a full boundary                         41.14
    integer                            :: IP, IPP    ! point counter                                    41.14
    integer                            :: IX         ! vertex number                                    41.14
    integer                            :: ISH        ! shift number                                     41.14
    !
    integer, dimension(:), allocatable :: blistot    ! list of all boundary vertices in ascending order
    integer, dimension(:), allocatable :: IARR1, IARR2 ! temporary array                                41.14
    !
    real                               :: d1         ! distance of a point to origin
    real                               :: d2         ! distance of another point to origin
    real                               :: xp, yp     ! coordinates of a boundary point                  41.14
    !
    character(80)                      :: msgstr     ! string to pass message
    character (len=8)                  :: PSNAME     ! name of output curve                             41.14
    !
    logical                            :: firstvert  ! indicate whether considered vertex is first vertex of boundary polygon
    logical                            :: found      ! indicates whether a new boundary part was found  41.14
    !
    type(celltype), dimension(:), pointer :: cell    ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face    ! datastructure for faces with their attributes
    type(verttype), dimension(:), pointer :: vert    ! datastructure for vertices with their attributes
    !
    type(OPSDAT), pointer :: OPSTMP                                   ! 41.14
    type XYPT                                                         ! 41.14
        real                :: X, Y                                   ! 41.14
        type(XYPT), pointer :: NEXTXY
    end type XYPT
    type(XYPT), target  :: FRST                                       ! 41.14
    type(XYPT), pointer :: CURR, TMP                                  ! 41.14
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanBpntlist')
    !
    ! if list of boundary vertices is already filled, return
    !
    if (allocated(blist)) return
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    vert(:)%atti(BINDX) = 0
    vert(:)%atti(BPOL)  = 0
    nbpt                = 0
    !
    ! determine total number of boundary vertices
    !
    nbptot = count(mask=vert(:)%atti(VMARKER)==1)
    !
    allocate(blistot(nbptot))
    blistot = 0
    !
    ! determine first boundary vertex nearest to the origin
    !
    kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VMARKER)==1)
    ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VMARKER)==1)
    !
    if ( kx(1) == ky(1) ) then
       !
       vc = kx(1)
       !
    else
       !
       d1 = sqrt((vert(kx(1))%attr(VERTX))**2+(vert(kx(1))%attr(VERTY))**2)
       d2 = sqrt((vert(ky(1))%attr(VERTX))**2+(vert(ky(1))%attr(VERTY))**2)
       !
       if ( d1 < d2 ) then
          vc = kx(1)
       else
          vc = ky(1)
       endif
       !
    endif
    !
    ! store first boundary vertex
    !
    vcf                 = vc
    blistot(1)          = vc
    nbpol               = 1
    firstvert           = .true.
    vert(vc)%atti(BPOL) = nbpol
    !
    nptemp = 0
    !
    ! store next subsequent boundary vertices in ascending order
    !
    k     = 1
    iface = 1
    !
    faceloop: do
       !
       if ( face(iface)%atti(FMARKER) == 1 ) then
          !
          if ( firstvert ) then
             !
             icell = face(iface)%atti(FACEC1)
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             ! pick up next vertex (counterclockwise counting of vertices is assumed)
             !
             vn = 0
             do j = 1, cell(icell)%nov
                if ( v(j) == vc ) then
                   vn = v(mod(j,cell(icell)%nov)+1)
                   exit
                endif
             enddo
             !
             if ( vn == 0 ) goto 10
             if ( vert(vn)%atti(VMARKER) /= 1 ) goto 10
             firstvert = .false.
             !
          endif
          !
          v1 = face(iface)%atti(FACEV1)
          v2 = face(iface)%atti(FACEV2)
          !
          if ( v1 == vc ) then
             !
             if ( v2 == vcf ) vc = vcf
             if ( any( v2 == blistot ) ) goto 10
             !
             k = k + 1
             blistot(k)          = v2
             vert(v2)%atti(BPOL) = nbpol
             vc = v2
             !
          elseif ( v2 == vc ) then
             !
             if ( v1 == vcf ) vc = vcf
             if ( any( v1 == blistot ) ) goto 10
             !
             k = k + 1
             blistot(k)          = v1
             vert(v1)%atti(BPOL) = nbpol
             vc = v1
             !
          elseif ( vc == vcf ) then        ! end of considered boundary polygon is found
             !
!ADC             exit faceloop
             if ( any( v1 == blistot ) .and. any( v2 == blistot ) ) goto 10
             !
             ! store number of boundary vertices for present polygon
             !
             nbpt(nbpol) = k - nptemp
             nptemp = k
             !
             ! take first vertex of next boundary polygon
             !
             vc                  = v1
             vcf                 = vc
             k = k + 1
             blistot(k)          = vc
             nbpol               = nbpol + 1
             firstvert           = .true.
             vert(vc)%atti(BPOL) = nbpol
             !
             ! give error if more than 10000 boundary polygons are found
             !
             if ( nbpol > 10000 ) call msgerr ( 2, ' More than 100 boundary polygons are found in grid' )
             !
          endif
          !
          if ( k == nbptot ) exit faceloop
          !
       endif
       !
 10    continue
       iface = iface + 1
       if ( iface > nfaces ) iface = 1
       !
    enddo faceloop
    !
    ! store number of boundary vertices for last polygon
    !
    nbpt(nbpol) = nbptot - nptemp
!ADC    nbpt(nbpol) = k - nptemp
    !
    ! check if list contains boundary vertices only
    !
    do j = 1, nbptot
!ADC    do j = 1, k
       !
       vc = blistot(j)
       if (vert(vc)%atti(VMARKER) /= 1) then
          write (msgstr, '(a,i4,a)') ' Vertex with index ',vc,' in boundary list is not a valid boundary point'
          call msgerr( 2, trim(msgstr) )
!ADC          vert(vc)%atti(VMARKER) = 1
       endif
       !
    enddo
    !
    ! determine maximum number of boundary vertices in set of polygons and allocate blist
    !
    maxnbp = maxval(nbpt)
    !
    if(.not.allocated(blist)) allocate (blist(maxnbp,nbpol), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanBpntlist: array blist ' )
       return
    endif
    blist = 0
    !
    ! fill blist in appropriate manner
    !
    k = 0
    !
    do j = 1, nbpol
       !
       do m = 1, nbpt(j)
          vc                   = blistot(k+m)
          blist(m,j)           = vc
          vert(vc)%atti(BINDX) = m
       enddo
       !
       k = k + nbpt(j)
       !
    enddo
    !
    deallocate(blistot)
    !
    !
    !  Add output curve corresponding to boundary                         41.14
    !
    ALLOCATE(OPSTMP)
    OPSTMP%PSTYPE = 'C'
    MIP = nbpt(1)
    OPSTMP%MIP = MIP
    ALLOCATE(OPSTMP%XP(MIP))
    ALLOCATE(OPSTMP%YP(MIP))
    OPSTMP%PSNAME = 'BOUNDARY'
    do m = 1, MIP
      vc = blist(m,1)
      OPSTMP%XP(m) = vert(vc)%attr(VERTX)
      OPSTMP%YP(m) = vert(vc)%attr(VERTY)
    enddo
    IF (ITEST.GE.10) WRITE (PRTEST, 104) 'BOUNDARY', MIP
104 format (' Generated output curve ', A8, ' with ', I4, ' vertices.')
!   ***** store number of points of the curve *****
    NULLIFY(OPSTMP%NEXTOPS)
    IF ( .NOT.LOPS ) THEN
       FOPS = OPSTMP
       COPS => FOPS
       LOPS = .TRUE.
    ELSE
       COPS%NEXTOPS => OPSTMP
       COPS => OPSTMP
    END IF
    !
    ! Determine highst value of VM
    !
    VMMAX = 0
    DO JBG = 1, nbpol
      DO IP = 1, nbpt(JBG)
        IX = blist(IP,JBG)
        VMMAX = MAX(VMMAX, vmark(IX))
      ENDDO
    ENDDO
!TEST    write (prtest, *) 'test VMMAX ', VMMAX, nbpol
    !
    ALLOCATE(IARR1(SUM(nbpt)))
    DO VM=1, VMMAX
      MIP = 0
      JJ = 0
      DO JBG = 1, nbpol
        !
        ! first boundary polygon is assumed an outer one
        ! (sea/mainland boundary) and hence, content of blist
        ! is ordered in counterclockwise manner
        !
        DO IP = 1, nbpt(JBG)
          IX = blist(IP,JBG)
          IF ( vmark(IX) == VM ) THEN
            MIP = MIP+1
            IARR1(MIP) = IP
            if (JJ==0) then
              JJ=JBG
            elseif (JJ/=JBG) then
              !call msgerr (1, 'Side is part of 2 boundaries')
            endif
          ENDIF
        ENDDO
      ENDDO
      !
      IF ( MIP/=0 ) THEN
        !
        ALLOCATE(IARR2(MIP))
        IARR2(1:MIP) = IARR1(1:MIP)
        ISH = 0
        DO IPP = 2, MIP
          IF ( IARR2(IPP)/=IARR2(IPP-1)+1 ) THEN
            ISH = IPP-1
            EXIT
          ENDIF
        ENDDO
        IARR2 = CSHIFT(IARR2,ISH)
!TEST        write (prtest, *) 'Shift ', ISH, MIP, IARR2(1)
        !
        ALLOCATE(OPSTMP)
        OPSTMP%PSTYPE = 'C'
        OPSTMP%MIP = MIP
        ALLOCATE(OPSTMP%XP(MIP))
        ALLOCATE(OPSTMP%YP(MIP))
        write (PSNAME, 101) VM
 101    format ('BOUND_',I2.2)
        OPSTMP%PSNAME = PSNAME
        DO IPP = 1, MIP
          IP = IARR2(IPP)
          IX = blist(IP,JJ)
          OPSTMP%XP(IPP) = vert(IX)%attr(VERTX)
          OPSTMP%YP(IPP) = vert(IX)%attr(VERTY)
        ENDDO
        DEALLOCATE(IARR2)
        IF (ITEST.GE.10) WRITE (PRTEST, 104) PSNAME, MIP
        NULLIFY(OPSTMP%NEXTOPS)
        IF ( .NOT.LOPS ) THEN
          FOPS = OPSTMP
          COPS => FOPS
          LOPS = .TRUE.
        ELSE
          COPS%NEXTOPS => OPSTMP
          COPS => OPSTMP
        END IF
        !
      ENDIF
    ENDDO
    DEALLOCATE(IARR1)
    !
end subroutine SwanBpntlist
