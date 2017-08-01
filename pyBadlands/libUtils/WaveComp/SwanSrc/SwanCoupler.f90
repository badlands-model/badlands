module swan_coupler_parallel

    use miscdata

    implicit none

    integer, save, allocatable :: CROSS(:)
    integer, save, allocatable :: BGRIDP(:)
    real, save,    allocatable :: BSPECS(:,:,:,:)
    real, save, allocatable :: AC1(:,:,:), COMPDA(:,:)
    real, save, allocatable    :: BLKND(:), BLKNDC(:), OURQT(:)

    ! SWAN model decomposition for SGFM framework
    integer :: Xsdiv, Ysdiv
    type REEF_WADecomp
        integer :: X_nb
        integer :: Y_nb
        integer :: minX_id
        integer :: minY_id
        integer :: maxX_id
        integer :: maxY_id
    end type REEF_WADecomp

end module swan_coupler_parallel
!=======================================================================
module swan_coupler_export

    implicit none

    real, dimension( :,:,: ), allocatable :: exportSWANfields

contains

    !=======================================================================
    subroutine export_swan( MIP, MXK, MYK, NVOQP, VOQ, VOQR, IONOD )

        use miscdata
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use M_BNDSPEC
        use OUTP_DATA
        use M_GENARR
        use M_PARALL

        integer :: MIP, MXK, MYK, IRQ, NVOQP
        integer :: VOQR(NMOVAR)
        integer :: IONOD(*)

        integer :: IV1, IV2, IV3, IV4, IV5, id, m, n, Isize, yik, xik, inerr, scenb

        real :: VOQ( MIP, NVOQP)
        real, dimension( MXC, MYC ) :: data

        if( NPROC == 1 )then
            Isize = MXCGL
        else
            if( MXCGL > MYCGL ) then
                Isize = MXC - IHALOX * IBLKAD( 1 )
            else
                Isize = MXC
            endif
        endif

        ! Average wave direction DIR
        IV1 = 13
        ! Significant wave height HSIGN
        IV2 = 10
        ! Average absolute wave period PER
        IV3 = 42
        ! Average wave lenght WLEN
        IV4 = 17
        ! RMS of orbital velocity amplitude at the bottom UBOT
        IV5 = 6

        m = 1
        n = 1
        call mpi_barrier( ocean_comm_world, inerr )

        do yik = 1, MYK
            id = ( yik - 1 )*MXK
            do xik = 1, MXK
                if( IONOD( id + xik) == INODE .or. NPROC == 1 )then
                    exportSWANfields( m, n, 1 ) = real( VOQ( id+ xik, VOQR(IV5)) )
                    exportSWANfields( m, n, 2 ) = real( VOQ( id+ xik, VOQR(IV1)) )
                    exportSWANfields( m, n, 3 ) = real( VOQ( id+ xik, VOQR(IV2)) )
                    exportSWANfields( m, n, 4 ) = real( VOQ( id+ xik, VOQR(IV3)) )
                    exportSWANfields( m, n, 5 ) = real( VOQ( id+ xik, VOQR(IV4)) )
                    m = m + 1
                    if( m > Isize )then
                        m = 1
                        n = n + 1
                    endif
                endif
            enddo
        enddo

        call mpi_barrier( ocean_comm_world, inerr )

        return

    end subroutine export_swan
    !=======================================================================

end module swan_coupler_export
!=======================================================================
module swan_coupler_functions

    use swan_coupler_parallel
    use swan_coupler_export
    use miscdata

    include 'mpif.h'

    integer :: nbX, nbY

    ! Imported spectral wave initial conditions
    type ImpSpecWa
        ! [hs] the significant wave height (in m).
        real :: hs
        ! [per] the characteristic period of the energy spectrum
        real :: per
        ! [dir] the peak wave direction
        real :: dir
    end type ImpSpecWa

    real, dimension(:,:), allocatable :: bathyfield, bathyHalo
    real, dimension(6) :: hcasts

contains

    !=======================================================================
    subroutine swan_decomposition_grid( nbPet, SWANdecomp )

        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use M_BNDSPEC
        use M_GENARR
        use M_PARALL

        implicit none
        integer :: ks, nbPet, inerr
        integer :: Isize, Jsize
        integer :: fdx, fdy, ldx, ldy
        integer :: minX, maxX, minY, maxY
        integer, dimension( nbPet ) :: minXT, maxXT, minYT, maxYT, nbXT, nbYT

        type( REEF_WADecomp ), dimension( nbPet ) :: SWANdecomp

        ldx = 0
        ldy = 0
        fdx = 0
        fdy = 0
        if( NPROC == 1 )then
            Isize = MXCGL
            Jsize = MYCGL
        else
            if( MXCGL > MYCGL ) then
                Isize = MXC - IHALOX * IBLKAD( 1 )
                Jsize = MYC
            else
                Isize = MXC
                Jsize = MYC - IHALOY * IBLKAD( 1 )
            endif
        endif

        if( MXF > 1 ) fdx = IHALOX
        if( MXL < MXCGL ) ldx = IHALOX
        if( MYF > 1 ) fdy = IHALOY
        if( MYL < MYCGL ) ldy = IHALOY


        nbX = Isize
        nbY = Jsize
        minX = MXF + fdx
        maxX = MXL - ldx
        minY = MYF + fdy
        maxY = MYL - ldy

        call mpi_allgather( nbX, 1, SWINT, nbXT, 1, SWINT, ocean_comm_world, inerr )
        call mpi_allgather( nbY, 1, SWINT, nbYT, 1, SWINT, ocean_comm_world, inerr )
        call mpi_allgather( minX, 1, SWINT, minXT, 1, SWINT, ocean_comm_world, inerr )
        call mpi_allgather( minY, 1, SWINT, minYT, 1, SWINT, ocean_comm_world, inerr )
        call mpi_allgather( maxX, 1, SWINT, maxXT, 1, SWINT, ocean_comm_world, inerr )
        call mpi_allgather( maxY, 1, SWINT, maxYT, 1, SWINT, ocean_comm_world, inerr )

        Xsdiv = 1
        Ysdiv = 1

        do ks = 1, NPROC
            SWANdecomp( ks )%X_nb = nbXT( ks )
            SWANdecomp( ks )%Y_nb = nbYT( ks )
            SWANdecomp( ks )%minX_id = minXT( ks )
            SWANdecomp( ks )%minY_id = minYT( ks )
            SWANdecomp( ks )%maxX_id = maxXT( ks )
            SWANdecomp( ks )%maxY_id = maxYT( ks )
            if( ks > 1 )then
                if( minXT( ks ) == maxXT( ks - 1 )+ 1 ) Xsdiv = Xsdiv + 1
                if( minYT( ks ) == maxYT( ks - 1 ) + 1 ) Ysdiv = Ysdiv + 1
            endif
        enddo
        call mpi_barrier( ocean_comm_world, inerr )
        allocate( bathyfield( nbX, nbY ) )
        allocate( bathyHalo( MXC, MYC ) )
        allocate( exportSWANfields( nbX, nbY, 5 ) )

    end subroutine swan_decomposition_grid
    !=======================================================================
    subroutine import_bathymetry

        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use M_BNDSPEC
        use M_GENARR
        use M_PARALL

        implicit none

        integer :: id, n, m, Isize, Jsize, Istr, Iend, Jstr, Jend, o, p, inerr
        integer :: lenght, l_bound, u_bound, left, right

        integer, dimension( mpi_status_size ) :: status
        real, dimension(:), allocatable :: sdl_1, sdl_2, sdl_3
        real, dimension(:), allocatable :: sdu_1, sdu_2, sdu_3
        real, dimension(:), allocatable :: rcl_1, rcl_2, rcl_3
        real, dimension(:), allocatable :: rcu_1, rcu_2, rcu_3

        ! Move values at 'present' time level 2 to 'old' time level 1.
        ! MCGRD = MXC*MYC+1-#masked cells.
        ! MXC = # cells x-dir in this tile including halox.
        ! MYC = # cells y-dir in this tile including haloy.
        ! COMPDA has only active wet points + 1.
        do id = 2, MCGRD
            compda( id, JDP1 ) = compda( id, JDP2 )
        enddo
        Isize = 0
        Jsize = 0

        if( nproc == 1 )then
            Istr=1
            Iend=MXCGL
            Jstr=1
            Jend=MYCGL
        else
            if( MXCGL > MYCGL )then
                if( INODE == 1 )then
                    Istr = 1
                else
                    Istr = IHALOX + 1
                endif
                Isize = MXC - IHALOX * IBLKAD(1)
                Iend = Istr + Isize - 1
                Jstr = 1
                Jend = MYC
            else
                if( INODE == 1 )then
                    Jstr = 1
                else
                    Jstr = IHALOY + 1
                endif
                Jsize = MYC - IHALOY * IBLKAD(1)
                Jend = Jstr + Jsize - 1
                Istr=1
                Iend=MXC
            endif
        endif
        ! We need to update halo points with new bathymetry arrays
        bathyHalo = 0.0
        p = 0
        do n = Jstr, Jend
            p = p + 1
            o = 0
            do m = Istr, Iend
                o = o + 1
                bathyHalo( m, n ) = bathyfield( o, p )
            enddo
        enddo

        ! Exchange halo cell data
        if( nproc > 1 )then

            if( MXCGL > MYCGL )then
                lenght = nbY
                l_bound = Istr
                u_bound = Iend
            else
                lenght = nbX
                l_bound = Jstr
                u_bound = Jend
            endif
            allocate( sdl_1( lenght ), sdl_2( lenght ), sdl_3( lenght ) )
            allocate( sdu_1( lenght ), sdu_2( lenght ), sdu_3( lenght ) )
            allocate( rcl_1( lenght ), rcl_2( lenght ), rcl_3( lenght ) )
            allocate( rcu_1( lenght ), rcu_2( lenght ), rcu_3( lenght ) )

            rcl_1 = 0.0_8
            rcl_2 = 0.0_8
            rcl_3 = 0.0_8
            rcu_1 = 0.0_8
            rcu_2 = 0.0_8
            rcu_3 = 0.0_8
            if( MXCGL > MYCGL )then
                sdl_1 = bathyHalo(l_bound,:)
                sdl_2 = bathyHalo(l_bound+1,:)
                sdl_3 = bathyHalo(l_bound+2,:)
                sdu_1 = bathyHalo(u_bound-2,:)
                sdu_2 = bathyHalo(u_bound-1,:)
                sdu_3 = bathyHalo(u_bound,:)
            else
                sdl_1 = bathyHalo(:,l_bound)
                sdl_2 = bathyHalo(:,l_bound+1)
                sdl_3 = bathyHalo(:,l_bound+2)
                sdu_1 = bathyHalo(:,u_bound-2)
                sdu_2 = bathyHalo(:,u_bound-1)
                sdu_3 = bathyHalo(:,u_bound)
            endif

            right = INODE
            left = INODE - 2
            if( right == nproc ) right = 0
            if( left < 0 ) left = nproc - 1

            call mpi_sendrecv( sdu_1, lenght, SWreal, right, 1232,  rcl_1, lenght, SWreal, left, 1232, &
                ocean_comm_world, status, inerr )
            call mpi_sendrecv( sdu_2, lenght, SWreal, right, 1234,  rcl_2, lenght, SWreal, left, 1234, &
                ocean_comm_world, status, inerr )
            call mpi_sendrecv( sdu_3, lenght, SWreal, right, 1235,  rcl_3, lenght, SWreal, left, 1235, &
                ocean_comm_world, status, inerr )

            call mpi_sendrecv( sdl_1, lenght, SWreal, left, 1231,  rcu_1, lenght, SWreal, right, 1231, &
                ocean_comm_world, status, inerr )
            call mpi_sendrecv( sdl_2, lenght, SWreal, left, 1236,  rcu_2, lenght, SWreal, right, 1236, &
                ocean_comm_world, status, inerr )
            call mpi_sendrecv( sdl_3, lenght, SWreal, left, 1237,  rcu_3, lenght, SWreal, right, 1237, &
                ocean_comm_world, status, inerr )

            ! Update bathymetry
            do n = 1, lenght
                if( MXCGL > MYCGL )then
                    if( INODE > 1 )then
                        bathyHalo(l_bound-3,n) = rcl_1(n)
                        bathyHalo(l_bound-2,n) = rcl_2(n)
                        bathyHalo(l_bound-1,n) = rcl_3(n)
                    endif
                    if( INODE < nproc )then
                        bathyHalo(u_bound+1,n) = rcu_1(n)
                        bathyHalo(u_bound+2,n) = rcu_2(n)
                        bathyHalo(u_bound+3,n) = rcu_3(n)
                    endif
                else
                    if( INODE > 1 )then
                        bathyHalo(n,l_bound-3) = rcl_1(n)
                        bathyHalo(n,l_bound-2) = rcl_2(n)
                        bathyHalo(n,l_bound-1) = rcl_3(n)
                    endif
                    if( INODE < nproc )then
                        bathyHalo(n,u_bound+1) = rcu_1(n)
                        bathyHalo(n,u_bound+2) = rcu_2(n)
                        bathyHalo(n,u_bound+3) = rcu_3(n)
                    endif
                endif
            enddo

            deallocate( sdl_1, sdl_2, sdl_3, sdu_1, sdu_2, sdu_3 )
            deallocate( rcl_1, rcl_2, rcl_3, rcu_1, rcu_2, rcu_3 )

        endif

        do n = 1, MYC
            do m = 1, MXC
                id = KGRPNT(m, n)
                if( id > 1 )then
                    compda( id, JDP2 ) = bathyHalo( m, n )
                endif
            enddo
        enddo

        return

    end subroutine import_bathymetry
    !=======================================================================
    subroutine swan_falsereadwind( COMPUT )

        use TIMECOMM
        use OCPCOMM1
        use OCPCOMM2
        use OCPCOMM3
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use OUTP_DATA
        use M_SNL4
        use M_GENARR
        use M_OBSTA
        use M_PARALL
        use SwanGriddata

        integer, parameter :: MXINCL = 10
        integer, parameter :: NVOTP  = 15
        integer, save     :: INCNUM(1:MXINCL) = 0
        integer, save     :: INCLEV = 1

        integer           :: IOSTAT = 0
        integer           :: ICNL4, ILAMBDA
        integer           :: ITMP1, ITMP2, ITMP3, ITMP4

        logical           :: MORE
        logical, save     :: LOBST = .false.

        integer        :: LREF, LREFDIFF, LRFRD
        real           :: POWN, DUM
        real           :: FD1, FD2, FD3, FD4
        real, allocatable :: RLAMBDA(:)

        character*6  QOVSNM
        character*40 QOVLNM
        real         QR(9)

        type(OBSTDAT), pointer :: OBSTMP
        type(OBSTDAT), save, pointer :: COBST

        type(OPSDAT), pointer :: OPSTMP

        type XYPT
            real                :: X, Y
            type(XYPT), pointer :: NEXTXY
        end type XYPT

        type(XYPT), target  :: FRST
        type(XYPT), pointer :: CURR, TMP

        type AUXT
            integer             :: I
            type(AUXT), pointer :: NEXTI
        end type AUXT
        type(AUXT), target  :: FRSTQ
        type(AUXT), pointer :: CURRQ, TMPQ

        type VEGPT
            integer              :: N
            real                 :: H, D, C
            type(VEGPT), pointer :: NEXTV
        end type VEGPT

        type(VEGPT), target  :: FRSTV
        type(VEGPT), pointer :: CURRV, TMPV

        logical :: found
        logical, save :: RUNMADE = .false.
        logical, save :: LOGCOM(1:6) = .false.

        integer   TIMARR(6)
        character PSNAME *8, PNAME *8, COMPUT *(*), PTYPE *1, DTTIWR *18
        integer, save :: IENT = 0       ! number of entries to this subr
        integer, allocatable :: IARR(:)
        integer IVOTP(NVOTP), INDX(1)
        data IVOTP /10, 11, 13, 15, 16, 17, 18, 19,          &
            28, 32, 33, 42, 43, 47, 48    /

        call swan_falseinitwind( AC2, SPCSIG, SPCDIR, KGRPNT )

        ! Read fake next compute
        COMPUT = 'COMP'
        RUNMADE = .true.
        NSTATM = 0
        TFINC = TINIC
        TIMCO = TINIC
        DT = 1.E10
        RDTIM = 0.
        NSTATC = 0
        MTC = 1
        NCOMPT = NCOMPT + 1
        if( NCOMPT > 50000 ) call MSGERR (2,   &
            'No more than 50000 COMPUTE commands are allowed')
        RCOMPT(NCOMPT,1) = REAL(NSTATC)
        RCOMPT(NCOMPT,2) = REAL(MTC)
        RCOMPT(NCOMPT,3) = TFINC
        RCOMPT(NCOMPT,4) = TINIC
        RCOMPT(NCOMPT,5) = DT
        ITERMX = MXITST

        return

    end subroutine swan_falsereadwind
    !=======================================================================
    subroutine swan_falseinitwind( AC2, SPCSIG, SPCDIR, KGRPNT )

        use OCPCOMM1
        use OCPCOMM2
        use OCPCOMM3
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use TIMECOMM
        use M_PARALL
        use SwanGriddata

        ! Initva parameters
        real SPCDIR(MDC,6), SPCSIG(MSC)
        logical :: STPNOW, EQCSTR
        integer    KGRPNT(MXC,MYC)
        integer JXMAX, JYMAX, NPTOT
        real       AC2(MDC,MSC,MCGRD)
        logical SINGLEHOT, PTNSUBGRD
        logical    KEYWIS, LERR
        character  RLINE *80
        integer IID, IUNITAC
        real    ACTMP(MDC), DIRTMP(MDC)
        logical EQREAL

        save       IENT
        data       IENT /0/

        U10 = hcasts( 4 )
        WDIP = hcasts( 5 )

        ! Convert (if necessary) WDIP from nautical degrees
        ! to cartesian degrees
        WDIP = DEGCNV (WDIP)

        if( IWIND == 0) IWIND = 4
        ALTMP = WDIP / 360.
        WDIP = PI2 * (ALTMP - NINT(ALTMP))

        return

    end subroutine swan_falseinitwind
    !=======================================================================
    subroutine swan_reinitialise

        use OCPCOMM2
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use TIMECOMM
        use OUTP_DATA
        use M_GENARR
        use M_PARALL
        use SwanGriddata
        use swan_coupler_parallel

        logical STPNOW
        integer ILEN, ISTAT, IF1, IL1
        character COMPUT *4
        character*20 NUMSTR, CHARS(1)
        character*80 MSGSTR

        call swan_initialisezero

        call swan_falsereadwind( COMPUT )

        ! Allocate some arrays meant for computation
        if( NUMOBS > 0 )then
            if( OPTG /= 5 )then
                ! Structured grid
                ILEN = 2*MCGRD
            else
                ! Unstructured grid
                ILEN = nfaces
            endif
            if(.not.allocated(CROSS)) allocate(CROSS(ILEN))
        else
            if (.not.allocated(CROSS)) allocate(CROSS(0))
        endif
        if (.not.allocated(BSPECS)) allocate(BSPECS(MDC,MSC,NBSPEC,2))
        if (.not.allocated(BGRIDP)) allocate(BGRIDP(6*NBGRPT))

        ! Do some preparations before computation
        call SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID, YCGRID, KGRPNT, &
            KGRBND, SPCDIR, SPCSIG )
        if (OPTG == 5) call SwanPrepComp ( CROSS )
        if (STPNOW()) return

        ! Check all possible flags and if necessary change
        call ERRCHK
        if (STPNOW()) return

        ! Initialisation of necessary grids for depth,
        !  current, wind and friction
        ISTAT = 0
!        ALOCMP = .true.
        if (ALOCMP.and.allocated(COMPDA)) deallocate(COMPDA)
        if (.not.allocated(COMPDA)) then
            allocate(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)
            ALOCMP = .false.
        endif
        if ( ISTAT/= 0 )then
            CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
            call TXPBLA(CHARS(1),IF1,IL1)
            MSGSTR =  'Allocation problem: array COMPDA and return code is '//   &
                CHARS(1)(IF1:IL1)
            call MSGERR ( 4, MSGSTR )
            return
        endif
        call SWRBC(COMPDA)

        ! Allocate AC1 in case of non-stationary situation or in case
        !  of using the S&L scheme

        if ( NSTATM == 1 .and. MXITNS > 1 .or. PROPSC == 3 ) then
            if (.not.allocated(AC1)) then
                allocate(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
            elseif(SIZE(AC1) == 0)then
                deallocate(AC1)
                allocate(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
            endif
            if ( ISTAT /= 0 ) then
                CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
                call TXPBLA(CHARS(1),IF1,IL1)
                MSGSTR = 'Allocation problem: array AC1 and return code is '//   &
                    CHARS(1)(IF1:IL1)
                call MSGERR ( 4, MSGSTR )
                return
            endif
            AC1 = 0.
        else
            if(.not.allocated(AC1)) allocate(AC1(0,0,0))
        endif

        if (LEVERR > MAXERR) then
            write (PRINTF, 6010) LEVERR
            if (LEVERR < 4) write (PRINTF, 6011)
6010        format(' ** No start of computation because of error level:',I3)
6011        format(' ** To ignore this error, change [maxerr] with the', ' SET command')

        else
            if (ITEST>=40) then
                if (NSTATC==1) then
                    write (PRINTF, '(" Type of computation: dynamic")')
                else
                    if (ONED) then
                        write (PRINTF, '(" Type of computation: static 1-D")')
                    else
                        write (PRINTF, '(" Type of computation: static 2-D")')
                    endif
                endif
            endif
        endif

    end subroutine swan_reinitialise
    !=======================================================================
    subroutine swan_reinitialise2(sideval)

        use OCPCOMM2
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use TIMECOMM
        use OUTP_DATA
        use M_GENARR
        use M_PARALL
        use SwanGriddata
        use SwanGridobjects
        use SwanCompdata
        use swan_coupler_parallel
        USE M_BNDSPEC

        logical STPNOW
        integer ILEN, ISTAT, IF1, IL1, sideval
        character COMPUT *4
        character*20 NUMSTR, CHARS(1)
        character*80 MSGSTR
        TYPE(BSDAT) , POINTER :: CURRBS
        real ANGLE

        ! Boundary setting
        if(sideval == 5)then
          ANGLE = 45.
        elseif(sideval == 7)then
          ANGLE = 135.
        elseif(sideval == 8)then
          ANGLE = -135.
        elseif(sideval == 6)then
          ANGLE = -45.
        elseif(sideval == 1)then
          ANGLE = 0.
        elseif(sideval == 3)then
          ANGLE = 90.
        elseif(sideval == 2)then
          ANGLE = 180.
        elseif(sideval == 4)then
          ANGLE = -90.
        endif

        call swan_initialisezero
        call SWBOUN2(ANGLE, hcasts, XCGRID, YCGRID, KGRPNT, XYTST, KGRBND)

        call swan_falsereadinput( COMPUT )

        !
        ! Allocate some arrays meant for computation
        if( NUMOBS > 0 )then
            if( OPTG /= 5 )then
                ! Structured grid
                ILEN = 2*MCGRD
            else
                ! Unstructured grid
                ILEN = nfaces
            endif
            if(.not.allocated(CROSS)) allocate(CROSS(ILEN))
        else
            if (.not.allocated(CROSS)) allocate(CROSS(0))
        endif
        if (.not.allocated(BSPECS)) allocate(BSPECS(MDC,MSC,NBSPEC,2))
        if (.not.allocated(BGRIDP)) allocate(BGRIDP(6*NBGRPT))

        ! CURRBS => FBS
        ! do
        !   NBS = CURRBS%NBS
        !   IF (NBS.EQ.-999) EXIT
        !   CURRBS%SPPARM(1:3) = hcasts(1:3) !SPPARM(1:4) = CURRBS%SPPARM(1:4)
        !   !CALL SSHAPE (BSPECS(1,1,NBS,1), SPCSIG, SPCDIR,  FSHAPE, DSHAPE)
        !   IF (.NOT.ASSOCIATED(CURRBS%NEXTBS)) EXIT
        !   CURRBS => CURRBS%NEXTBS
        ! enddo

        ! Do some preparations before computation
        call SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID, YCGRID, KGRPNT, &
            KGRBND, SPCDIR, SPCSIG )
        if (OPTG == 5) call SwanPrepComp ( CROSS )

        if (STPNOW()) return

        ! Check all possible flags and if necessary change
        call ERRCHK
        if (STPNOW()) return

        ! Initialisation of necessary grids for depth,
        !  current, wind and friction
        ISTAT = 0
        ALOCMP = .true.
        if (ALOCMP.and.allocated(COMPDA)) deallocate(COMPDA)
        if (.not.allocated(COMPDA)) then
            allocate(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)
            ALOCMP = .false.
        endif
        if ( ISTAT/= 0 )then
            CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
            call TXPBLA(CHARS(1),IF1,IL1)
            MSGSTR =  'Allocation problem: array COMPDA and return code is '//   &
                CHARS(1)(IF1:IL1)
            call MSGERR ( 4, MSGSTR )
            return
        endif
        call SWRBC(COMPDA)

        ! Allocate AC1 in case of non-stationary situation or in case
        ! of using the S&L scheme

        if ( NSTATM == 1 .and. MXITNS > 1 .or. PROPSC == 3 ) then
            if (.not.allocated(AC1)) then
                allocate(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
            elseif(SIZE(AC1) == 0)then
                deallocate(AC1)
                allocate(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
            endif
            if ( ISTAT /= 0 ) then
                CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
                call TXPBLA(CHARS(1),IF1,IL1)
                MSGSTR = 'Allocation problem: array AC1 and return code is '//   &
                    CHARS(1)(IF1:IL1)
                call MSGERR ( 4, MSGSTR )
                return
            endif
            AC1 = 0.
        else
            if(.not.allocated(AC1)) allocate(AC1(0,0,0))
        endif

        CALL SWSYNC

        if (LEVERR > MAXERR) then
            write (PRINTF, 6010) LEVERR
            if (LEVERR < 4) write (PRINTF, 6011)
6010        format(' ** No start of computation because of error level:',I3)
6011        format(' ** To ignore this error, change [maxerr] with the', ' SET command')

        else
            if (ITEST>=40) then
                if (NSTATC==1) then
                    write (PRINTF, '(" Type of computation: dynamic")')
                else
                    if (ONED) then
                        write (PRINTF, '(" Type of computation: static 1-D")')
                    else
                        write (PRINTF, '(" Type of computation: static 2-D")')
                    endif
                endif
            endif
        endif

    end subroutine swan_reinitialise2
    !=======================================================================
    subroutine swan_falsereadinput( COMPUT )

        use TIMECOMM
        use OCPCOMM1
        use OCPCOMM2
        use OCPCOMM3
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use OUTP_DATA
        use M_SNL4
        use M_GENARR
        use M_OBSTA
        use M_PARALL
        use SwanGriddata

        integer, parameter :: MXINCL = 10
        integer, parameter :: NVOTP  = 15
        integer, save     :: INCNUM(1:MXINCL) = 0
        integer, save     :: INCLEV = 1

        integer           :: IOSTAT = 0
        integer           :: ICNL4, ILAMBDA
        integer           :: ITMP1, ITMP2, ITMP3, ITMP4

        logical           :: MORE
        logical, save     :: LOBST = .false.

        integer        :: LREF, LREFDIFF, LRFRD
        real           :: POWN, DUM
        real           :: FD1, FD2, FD3, FD4
        real, allocatable :: RLAMBDA(:)

        character*6  QOVSNM
        character*40 QOVLNM
        real         QR(9)

        type(OBSTDAT), pointer :: OBSTMP
        type(OBSTDAT), save, pointer :: COBST

        type(OPSDAT), pointer :: OPSTMP

        type XYPT
            real                :: X, Y
            type(XYPT), pointer :: NEXTXY
        end type XYPT

        type(XYPT), target  :: FRST
        type(XYPT), pointer :: CURR, TMP

        type AUXT
            integer             :: I
            type(AUXT), pointer :: NEXTI
        end type AUXT
        type(AUXT), target  :: FRSTQ
        type(AUXT), pointer :: CURRQ, TMPQ

        type VEGPT
            integer              :: N
            real                 :: H, D, C
            type(VEGPT), pointer :: NEXTV
        end type VEGPT

        type(VEGPT), target  :: FRSTV
        type(VEGPT), pointer :: CURRV, TMPV

        logical :: found
        logical, save :: RUNMADE = .false.
        logical, save :: LOGCOM(1:6) = .false.

        integer   TIMARR(6)
        character PSNAME *8, PNAME *8, COMPUT *(*), PTYPE *1, DTTIWR *18
        integer, save :: IENT = 0       ! number of entries to this subr
        integer, allocatable :: IARR(:)
        integer IVOTP(NVOTP), INDX(1)
        data IVOTP /10, 11, 13, 15, 16, 17, 18, 19,          &
            28, 32, 33, 42, 43, 47, 48    /

        call swan_falseinitinput( AC2, SPCSIG, SPCDIR, KGRPNT )

        ! Read fake next compute
        COMPUT = 'COMP'
        RUNMADE = .true.
        NSTATM = 0
        TFINC = TINIC
        TIMCO = TINIC
        DT = 1.E10
        RDTIM = 0.
        NSTATC = 0
        MTC = 1
        NCOMPT = NCOMPT + 1
        if( NCOMPT > 50000 ) call MSGERR (2,   &
            'No more than 50000 COMPUTE commands are allowed')
        RCOMPT(NCOMPT,1) = REAL(NSTATC)
        RCOMPT(NCOMPT,2) = REAL(MTC)
        RCOMPT(NCOMPT,3) = TFINC
        RCOMPT(NCOMPT,4) = TINIC
        RCOMPT(NCOMPT,5) = DT
        ITERMX = MXITST

        return

    end subroutine swan_falsereadinput
    !=======================================================================
    subroutine swan_falseinitinput( AC2, SPCSIG, SPCDIR, KGRPNT )

        use OCPCOMM1
        use OCPCOMM2
        use OCPCOMM3
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use TIMECOMM
        use M_PARALL
        use SwanGriddata

        ! Initva parameters
        real SPCDIR(MDC,6), SPCSIG(MSC)
        logical :: STPNOW, EQCSTR
        integer    KGRPNT(MXC,MYC)
        integer JXMAX, JYMAX, NPTOT
        real       AC2(MDC,MSC,MCGRD)
        logical SINGLEHOT, PTNSUBGRD
        logical    KEYWIS, LERR
        character  RLINE *80
        integer IID, IUNITAC
        real    ACTMP(MDC), DIRTMP(MDC)
        logical EQREAL

        save       IENT
        data       IENT /0/

        ! Read fake new initial values for borders
!         call IGNORE ('COND')
!         LERR =  .false.
!         if (MXC <= 0 .and. OPTG /= 5)then
!             call MSGERR (2, 'command fake INIT should follow CGRID')
!             LERR = .true.
!         endif
!         if (MCGRD <= 1 .and. nverts <= 0) then
!             call MSGERR (2,'command fake INIT should follow READ BOT or READ UNSTRUC')
!             LERR = .true.
!         endif
        ICOND = 2
!         ! Hsign
        ! SPPARM(1)= hcasts( 1 )
!         if (KEYWIS('MEAN')) then
!             if (FSHAPE > 0) FSHAPE=-FSHAPE
!         elseif (KEYWIS('PEAK')) then
!             if (FSHAPE < 0) FSHAPE=-FSHAPE
!         endif
!         ! Period
        ! SPPARM(2) = hcasts( 2 )
!         if (SPPARM(2) <= 0.) call MSGERR (2, 'Period must be >0')
!         if (SPPARM(2) > PI2/SPCSIG(1)) then
!             call MSGERR (2,'Inc. freq. lower than lowest spectral freq.')
!             write(PRINTF,11) 1./SPPARM(2),' < ',SPCSIG(1)/PI2,' lowest'
!         endif
!         if (SPPARM(2) < PI2/SPCSIG(MSC)) then
!             call MSGERR (2,'Inc. freq. higher than highest spectral freq.')
!             write(PRINTF,11) 1./SPPARM(2),' > ',SPCSIG(MSC)/PI2,'highest'
!         endif
! 11      format(' Inc. freq. = ',F9.5, A3, F9.5, '=',A7,' freq')
!         ! Dir
        ! SPPARM( 3 ) = hcasts( 3 )

        ! Give boundaries of region where the initial condition applies
        ! IX1 = 0
        ! IX2 = MXC-1
        ! IY1 = 0
        ! IY2 = MYC-1
        ! if (.not.LERR)then
        !     call SSHAPE (AC2(1,1,1), SPCSIG, SPCDIR, FSHAPE, DSHAPE)
        !     do IX = IX1+1, IX2+1
        !         do IY = IY1+1, IY2+1
        !             INDX = KGRPNT(IX,IY)
        !             if (INDX.GT.1)call SINTRP (1., 0., AC2(1,1,1), AC2(1,1,1), &
        !                 AC2(1,1,INDX),SPCDIR,SPCSIG)
        !         enddo
        !     enddo
        !     ! Reset action density AC2(*,*,1) to 0
        !     do ID = 1, MDC
        !         do IS = 1, MSC
        !             AC2(ID,IS,1) = 0.
        !         enddo
        !     enddo
        ! endif

        U10 = hcasts( 4 )
        WDIP = hcasts( 5 )

        ! Convert (if necessary) WDIP from nautical degrees
        ! to cartesian degrees
        WDIP = DEGCNV (WDIP)

        if( IWIND == 0) IWIND = 4
        ALTMP = WDIP / 360.
        WDIP = PI2 * (ALTMP - NINT(ALTMP))

        return

    end subroutine swan_falseinitinput
    !=======================================================================
    subroutine swan_initialisezero !( AC2, SPCSIG, SPCDIR, KGRPNT )

        use OCPCOMM1
        use OCPCOMM2
        use OCPCOMM3
        use OCPCOMM4
        use SWCOMM1
        use SWCOMM2
        use SWCOMM3
        use SWCOMM4
        use TIMECOMM
        use M_GENARR
        use OUTP_DATA
        use M_OBSTA
        use M_SNL4
        use M_PARALL
        use SwanGriddata

        integer :: INDX, ID, IS

        do INDX = 1, MCGRD
            do ID = 1, MDC
                do IS = 1, MSC
                    AC2(ID,IS,INDX) = 0.
                enddo
            enddo
        enddo
        ICOND = 2

        return

    end subroutine swan_initialisezero
    !=======================================================================

end module swan_coupler_functions
