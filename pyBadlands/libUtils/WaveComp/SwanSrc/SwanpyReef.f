!
!     SWAN main program and miscellaneous routines
!
!     Contents of this file:
!
!     SWAN:   Main program
!     SWMAIN: Calling SWINIT, SWREAD, SWCOMP and SWOUTP
!     SWINIT: Initialize several variables and arrays
!     SWPREP: Do some preparations before computation is started
!     SPRCON: Execution of some tests on the given model description
!     SWRBC
!     SVALQI
!     SINUPT
!     SINBTG
!     SINCMP
!     WRTEST
!     ERRCHK
!     SNEXTI
!     RBFILE: Read boundary spectra from one file                         40.00
!     RESPEC: Read one 1-d OR 2-d boundary spectrum from file, and        40.00
!             transform to internal SWAN spectral resolution              40.00
!     FLFILE: Update boundary conditions, update nonstationary input      40.00
!             fields                                                      40.00
!     SWINCO
!     SWCLME: Clean memory
!
!***********************************************************************

!
! IN CASE OF COUPLING WITH OCEAN MODULE
! CHANGE AT LINE 2502 THE 0 to 1
!              UU  = SVALQI (XP, YP, 2, UXB, 1 ,IX ,IY)
!              VV  = SVALQI (XP, YP, 3, UYB, 1 ,IX ,IY)
!
      SUBROUTINE SWAN_INITIALIZE ( MyComm, inputname, infoname )

      USE TIMECOMM                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.51
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
      USE swan_coupler_functions
      use swan_coupler_parallel

      IMPLICIT NONE


      CHARACTER*80, intent(in) :: inputname, infoname
      INTEGER, INTENT(in) :: MyCOMM
      INTEGER   IUNIT                                                     34.01
      INTEGER   IOSTAT, IT0, IT, SAVITE, ILEN                             40.30
      INTEGER   INERR                                                     40.31
      INTEGER   IERR                                                      40.31 40.00
      INTEGER   ISTAT, IF1, IL1                                           40.41
      CHARACTER PTYPE, PNAME *8, COMPUT *4, DTTIWR*18                     40.00
      CHARACTER*20 NUMSTR, CHARS(1)                                       40.41
      CHARACTER*80 MSGSTR                                                 40.41
      LOGICAL   LOPEN
      LOGICAL STPNOW                                                      34.01

!     --- initialize MPI
      CALL SWINITMPI( MyComm )

!     --- initialize various data
      LEVERR=0                                                            40.23
      MAXERR=1                                                            34.01
      ITRACE=0                                                            40.23
      INERR =0                                                            40.31

      CALL SWINIT (INERR, inputname, infoname)                                        40.31 34.01
      IF (INERR.GT.0) RETURN                                              34.01
      IF (STPNOW()) RETURN                                                34.01
      COMPUT = '    '

!       --- read and process user commands
      CALL SWREAD (COMPUT)                                              40.31 30.90
      IF (STPNOW()) RETURN                                              34.01

      IF (COMPUT.EQ.'STOP') THEN                                        40.13
          IUNIT  = 0                                                      40.13
          IOSTAT = 0                                                      40.13
          FILENM = 'norm_end'                                             40.13
          IF  (INODE.EQ.MASTER) THEN
            CALL FOR (IUNIT, FILENM, 'UF', IOSTAT)                          40.13
            WRITE (IUNIT, *) ' Normal end of run ', PROJNR                  40.13
          ENDIF
          GOTO 900                                                        40.13
      ENDIF

!       --- allocate some arrays meant for computation                    40.31
        IF (NUMOBS .GT. 0) THEN
           IF (OPTG.NE.5) THEN                                            40.80
!             structured grid                                             40.80
              ILEN = 2*MCGRD                                              40.80
           ELSE                                                           40.80
!             unstructured grid                                           40.80
              ILEN = nfaces                                               40.80
           ENDIF                                                          40.80
           IF (.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(ILEN))               40.80 40.31
        ELSE
           IF (.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(0))                  40.80 40.31
        ENDIF                                                             34.01
        IF (.NOT.ALLOCATED(BSPECS)) ALLOCATE(BSPECS(MDC,MSC,NBSPEC,2))    40.31
        IF (.NOT.ALLOCATED(BGRIDP)) ALLOCATE(BGRIDP(6*NBGRPT))            40.31

!       --- do some preparations before computation                       40.31
        CALL SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID, YCGRID, KGRPNT,     40.31
     &                KGRBND, SPCDIR, SPCSIG )                            40.31
        IF (OPTG.EQ.5) CALL SwanPrepComp ( CROSS )                        40.80
        IF (STPNOW()) RETURN                                              40.80

!       --- check all possible flags and if necessary change
!           if option is not correct
        CALL ERRCHK                                                       30.60
        IF (STPNOW()) RETURN

!       --- initialisation of necessary grids for depth,
!           current, wind and friction
        IF (ALOCMP.AND.ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)              40.97
        IF (.NOT.ALLOCATED(COMPDA)) THEN                                  40.97
           ALLOCATE(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)                      40.97 40.41 40.31
           ALOCMP = .FALSE.                                               40.97
        END IF                                                            40.97
        IF ( ISTAT.NE.0 ) THEN                                            40.41
           CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                           40.41
           CALL TXPBLA(CHARS(1),IF1,IL1)                                  40.41
           MSGSTR =                                                       40.41
     &         'Allocation problem: array COMPDA and return code is '//   40.41
     &         CHARS(1)(IF1:IL1)                                          40.41
           CALL MSGERR ( 4, MSGSTR )                                      40.41
           RETURN                                                         40.41
        END IF                                                            40.41

        CALL SWRBC(COMPDA)                                                40.31

!       --- allocate AC1 in case of non-stationary situation or in case   40.31
!           of using the S&L scheme                                       40.31
        IF ( NSTATM.EQ.1 .AND. MXITNS.GT.1 .OR. PROPSC.EQ.3 ) THEN        40.31
           IF (.NOT.ALLOCATED(AC1)) THEN                                  40.41 40.31
              ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)                     40.41
           ELSE IF (SIZE(AC1).EQ.0) THEN                                  40.41
              DEALLOCATE(AC1)                                             40.41
              ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)                     40.41
           END IF                                                         40.41
           IF ( ISTAT.NE.0 ) THEN                                         40.41
              CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                        40.41
              CALL TXPBLA(CHARS(1),IF1,IL1)                               40.41
              MSGSTR =                                                    40.41
     &            'Allocation problem: array AC1 and return code is '//   40.41
     &            CHARS(1)(IF1:IL1)                                       40.41
              CALL MSGERR ( 4, MSGSTR )                                   40.41
              RETURN                                                      40.41
           END IF                                                         40.41
           AC1 = 0.                                                       40.31
        ELSE                                                              40.31
           IF(.NOT.ALLOCATED(AC1)) ALLOCATE(AC1(0,0,0))                   40.31
        ENDIF

        IF (LEVERR.GT.MAXERR) THEN                                        40.00

          WRITE (PRINTF, 6010) LEVERR
          IF (LEVERR.LT.4) WRITE (PRINTF, 6011)                           30.72
 6010     FORMAT(' ** No start of computation because of error level:'
     &      ,I3)
 6011     FORMAT(' ** To ignore this error, change [maxerr] with the',    30.72
     &           ' SET command')                                          30.72
        ENDIF

          IF (ITEST.GE.40) THEN                                           40.00
            IF (NSTATC.EQ.1) THEN                                         33.08
              WRITE (PRINTF, '(" Type of computation: dynamic")')         32.02
            ELSE                                                          32.02
              IF (ONED) THEN                                              32.02
                WRITE (PRINTF, '(" Type of computation: static 1-D")')    32.02
              ELSE                                                        32.02
                WRITE (PRINTF, '(" Type of computation: static 2-D")')    32.02
              ENDIF                                                       32.02
            ENDIF                                                         32.02
          ENDIF

          IF (NSTATC.EQ.1) THEN                                           40.00
            IT0 = 0                                                       40.00
            IF (ICOND.EQ.1) THEN                                          40.00
!             --- compute default initial conditions
              CALL SWINCO ( AC2   , COMPDA, XCGRID, YCGRID,               40.31
     &                      KGRPNT, SPCDIR, SPCSIG, XYTST )               40.31

!             --- reset ICOND to prevent second computation of
!                 initial condition
              ICOND = 0                                                   40.00
            ENDIF
          ELSE
            IT0 = 1
          ENDIF

 900  CONTINUE

      RETURN
      END SUBROUTINE SWAN_INITIALIZE
!***********************************************************************
!                                                                      *
      SUBROUTINE SWAN_RUN ( hindcast )
!                                                                      *
!***********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.51
      USE M_GENARR                                                        40.31
      USE M_BNDSPEC
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
      USE swan_coupler_functions
      use swan_coupler_parallel

      IMPLICIT NONE
!
      REAL,  intent(in) :: hindcast( 6 )
      INTEGER   IUNIT
      INTEGER   IOSTAT, IT0, IT, SAVITE, ILEN
      INTEGER   INERR
      INTEGER   IERR
      INTEGER   ISTAT, IF1, IL1
      CHARACTER PTYPE, PNAME *8, COMPUT *4, DTTIWR*18
      CHARACTER*20 NUMSTR, CHARS(1)
      CHARACTER*80 MSGSTR
      LOGICAL   LOPEN
      LOGICAL STPNOW

      TYPE(BSDAT) , POINTER :: CURRBS                                     40.31

          hcasts = hindcast

          LEVERR = 0
!         --- synchronize nodes                                           40.30
          CALL SWSYNC                                                     40.30
          IF (STPNOW()) RETURN                                            40.30

C          OPEN(UNIT=19485, FILE='ac2.txt', STATUS='NEW')
C          WRITE(UNIT=19485,FMT=*) AC2
C          CLOSE(UNIT=19485)
C          OPEN(UNIT=19487, FILE='BSPECS.txt', STATUS='NEW')
C          WRITE(UNIT=19487,FMT=*) BSPECS
C          CLOSE(UNIT=19487)
C          OPEN(UNIT=19488, FILE='COMPDA.txt', STATUS='NEW')
C          WRITE(UNIT=19488,FMT=*) COMPDA(1,JHSIBC)
C          CLOSE(UNIT=19488)

          IT0 = 1
!          --- loop over time steps                                        40.00
          DO 500 IT = IT0, MTC                                            40.00

            IF (LEVERR.GT.MAXERR) THEN                                    30.82
              WRITE (PRINTF, 6030) LEVERR                                 30.82
              IF (LEVERR.LT.4) WRITE (PRINTF, 6011)                       40.30 30.82
 6011     FORMAT(' ** To ignore this error, change [maxerr] with the',
     &           ' SET command')
 6030         FORMAT(' ** No continuation of computation because ',       30.82
     &               'of error level:',I3)                                30.82
              EXIT                                                        40.30
            ENDIF                                                         30.82

!           --- synchronize nodes
            CALL SWSYNC                                                   40.30
            IF (STPNOW()) RETURN                                          40.30

!           --- update boundary conditions and input fields
            CALL SNEXTI ( BSPECS, BGRIDP, COMPDA, AC1   , AC2   ,         40.31
     &                    SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT,         40.31
     &                    XYTST , DEPTH , WLEVL , FRIC  , UXB   ,         40.31
     &                    UYB   , NPLAF , WXI   , WYI   )                 40.55 40.31
            IF (STPNOW()) RETURN                                          34.01

!           --- synchronize nodes
            CALL SWSYNC                                                   40.30

            IF (STPNOW()) RETURN                                          40.30
            IF (COMPUT.NE.'NOCO' .AND. IT.GT.0) THEN                      40.00
              SAVITE = ITEST                                              30.21
              IF (ICOTES .GT. ITEST) ITEST = ICOTES

!             --- compute action density for current time step
              IF (OPTG.NE.5) THEN                                         40.80
!                structured grid                                          40.80
                 CALL SWCOMP( AC1   , AC2   , COMPDA, SPCDIR, SPCSIG,     40.31
     &                        XYTST , IT    , KGRPNT, XCGRID, YCGRID,     40.31
     &                        CROSS )                                     40.31
              ELSE                                                        40.80
!                unstructured grid                                        40.80
                 CALL SwanCompUnstruc ( AC2   , AC1   , COMPDA,           40.80
     &                                  SPCSIG, SPCDIR, XYTST ,           40.80
     &                                  CROSS , IT    )                   40.80
              ENDIF                                                       40.80

              IF (STPNOW()) RETURN                                        34.01
!             --- set ICOND=4 for stationary computation, for next
!                 (stationary) COMPUTE command                            40.13
              ICOND = 4                                                   40.13

!             --- check whether computed significant wave height at       32.01
!                 boundary differs from prescribed value given in         32.01
!                 boundary command values of incident Hs                  32.01
              IF ( BNDCHK ) THEN                                          32.01
                CALL HSOBND ( AC2, SPCSIG, COMPDA(1,JHSIBC), KGRPNT )     40.31
              ENDIF                                                       32.01
              ITEST = SAVITE                                              30.21
            ENDIF

            IF ( IT.EQ.IT0 .AND. .NOT.ALLOCATED(OURQT) ) THEN             40.51 40.31
               ALLOCATE (OURQT(MAX_OUTP_REQ))                             40.51 40.30
               OURQT = -9999.                                             40.51 40.30
            ENDIF                                                         40.00

            SAVITE = ITEST                                                30.21
            IF (IOUTES .GT. ITEST) ITEST = IOUTES

!           --- synchronize nodes
            CALL SWSYNC                                                   40.30
            IF (STPNOW()) RETURN                                          40.30

!           --- couple models during output computations
!            IF (IT.gt.0) THEN
              CALL SWOUTP ( AC2, SPCSIG,SPCDIR,COMPDA,XYTST ,
     &                      KGRPNT, XCGRID, YCGRID, KGRBND, OURQT )

!            END IF
            IF (STPNOW()) RETURN                                          40.30
!
            IF (ERRPTS.GT.0) REWIND(ERRPTS)                               30.50
            ITEST = SAVITE                                                30.21

!           --- update time
            IF (NSTATM.EQ.1) THEN                                         40.00
              IF (NSTATC.EQ.1 .AND. IT.LT.MTC) TIMCO = TIMCO + DT         40.00
              CHTIME = DTTIWR(ITMOPT, TIMCO)                              40.00
              IF (NSTATC.EQ.1) WRITE (PRINTF, 222) CHTIME, TIMCO          40.00
 222          FORMAT(' Time of computation ->  ',A,' in sec:', F12.0)     40.00
            ENDIF                                                         40.00

 500      CONTINUE

          call swan_reinitialise2(int(hcasts(6)))

      RETURN
      END SUBROUTINE SWAN_RUN
!***********************************************************************
!                                                                      *
      SUBROUTINE SWAN_FINALIZE
!                                                                      *
!***********************************************************************
      USE TIMECOMM                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.51
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
      USE swan_coupler_functions
      use swan_coupler_parallel
      use swan_coupler_functions

!
!     Local variables
!
!     BLKND :     array giving node number per subdomain
!     BLKNDC:     auxiliary array for collecting the array BLKND
!
      LOGICAL STPNOW
      LOGICAL LOPEN
      INTEGER MyStatus

      DO IUNIT=1,HIOPEN                                                   34.01
        INQUIRE(UNIT=IUNIT,OPENED=LOPEN)                                  34.01
        IF (LOPEN.AND.IUNIT.NE.PRINTF) CLOSE(IUNIT)                       40.30
      END DO                                                              34.01

!     --- collect contents of individual process files for                40.30
!         output requests in case of parallel computation                 40.30

      CALL SWSYNC                                                         40.30
      IF ( PARLL ) THEN                                                   40.30
         ALLOCATE (BLKND(MXC*MYC))                                        40.51 40.41
         BLKND = REAL(INODE)                                              40.41
         IF ( IAMMASTER ) THEN                                            40.95
            ALLOCATE(BLKNDC(MXCGL*MYCGL))                                 40.51 40.41
            BLKNDC = 0.                                                   40.96
         END IF
         CALL SWCOLLECT ( BLKNDC, BLKND, .TRUE. )                         40.51 40.41
         IF (STPNOW()) RETURN                                             40.41
         IF ( IAMMASTER ) THEN                                            40.95 40.30
            CALL SWCOLOUT ( OURQT, BLKNDC )                               40.51 40.41 40.30
            DEALLOCATE(BLKNDC)                                            40.41 40.30
         END IF                                                           40.30
      END IF                                                              40.30

      INQUIRE(UNIT=PRINTF,OPENED=LOPEN)                                   40.30
      IF (LOPEN) CLOSE(PRINTF)                                            40.30

!     --- deallocate all allocated arrays                                 40.31
      IF (ALLOCATED(AC1   )) DEALLOCATE(AC1   )                           40.31
      IF (ALLOCATED(BGRIDP)) DEALLOCATE(BGRIDP)                           40.31
      IF (ALLOCATED(BSPECS)) DEALLOCATE(BSPECS)                           40.31
      IF (ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)                           40.31
      IF (ALLOCATED(CROSS )) DEALLOCATE(CROSS )                           40.31
      IF (ALLOCATED(OURQT )) DEALLOCATE(OURQT )                           40.51 40.30
      IF (ALLOCATED(BLKND )) DEALLOCATE(BLKND )                           40.41

      CALL SWCLME                                                         40.31
C
C      CALL SWEXITMPI                                                      40.30

! THIS IS WHERE WE INITIALIZE ESFM COUPLING

      RETURN

      END SUBROUTINE SWAN_FINALIZE
!***********************************************************************
!                                                                      *
      SUBROUTINE SWINIT (INERR, inputname, infoname )                                           40.31 34.01
!                                                                      *
!***********************************************************************
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE TIMECOMM                                                        40.41
      USE M_SNL4                                                          40.17
      USE M_BNDSPEC                                                       40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80

      INTEGER INERR, IGRID, MXOUTAR

      LOGICAL STPNOW                                                      34.01
      character :: inputname*40, infoname*40

      VERTXT = BLANK                                                      40.03
      VERNUM = 40.85
      WRITE (VERTXT, '(F5.2)') VERNUM                                     40.03
!      CALL BUGFIX ('A')
!
      CALL OCPINI ('swaninit', inputname, infoname, .TRUE.,INERR)                              34.01
      IF (INERR.GT.0) RETURN                                              34.01
      IF (STPNOW()) RETURN                                                34.01
!
      WRITE (PRINTF, 6010) VERTXT                                         40.03
 6010 FORMAT (/,20X,'---------------------------------------',
     &        /,20X,'                 SWAN',
     &        /,20X,'SIMULATION OF WAVES IN NEAR SHORE AREAS',
     &        /,20X,'         VERSION NUMBER ', A,                        40.03
     &        /,20X,'---------------------------------------',//)
!
!      IF (SCREEN.NE.PRINTF.AND.IAMMASTER) WRITE (SCREEN,6020)             40.95 40.30
 6020 FORMAT (/, ' SWAN is preparing computation',/)
!
!     ***** initial values for common variables *****
!     ***** names *****
      PROJID = 'SWAN'
      PROJNR = BLANK
      PROJT1 = BLANK
      PROJT2 = BLANK
      PROJT3 = BLANK
      FNEST  = BLANK
      FBCR   = BLANK
      FBCL   = BLANK
      UH     = 'm'
      UV     = 'm/s'
      UT     = 'sec'
      UL     = 'm'
      UET    = 'm3/s'
      UDI    = 'degr'
      UST    = 'm2/s2'
      UF     = 'N/m2'
      UP     = 'W/m'
      UAP    = 'W/m2'
      UDL    = 'm2/s'
!     ***** physical parameters *****
      GRAV   = 9.81
      WLEV   = 0.
      CASTD  = 0.               ! const. air-sea temp diff                40.03
      CDCAP  = 99999.
      PI     = 4.*ATAN(1.)                                                40.31
      PI2    = 2.*PI
      UNDFLW = 1.E-15
      DNORTH = 90.                                                        30.72
      DEGRAD = PI/180.
      RHO    = 1025.
!     power of tail in spectrum, 1: E with f, 2: E with k,
!                                3: A with f, 4: A with k
      PWTAIL(1) = 4.
      PWTAIL(2) = 2.5
      PWTAIL(3) = PWTAIL(1)+1.                                            30.72
      PWTAIL(4) = 3.
!     ***** number of computational grid points ****                      30.60
      MCGRD   = 1                                                         30.60
      MCGRDGL = 1                                                         40.30
      NGRBND  = 0                                                         40.00
      NGRBGL  = 0                                                         40.30
      nverts  = 0                                                         40.80
      nvertsg = 0                                                         40.95
!     time of computation                                                 40.00
      TIMCO = -1.E10                                                      40.00
      CHTIME = '    '                                                     40.00
!     boundary conditions                                                 40.00
      NBFILS = 0                                                          40.00
      NBSPEC = 0                                                          40.00
      NBGRPT = 0                                                          40.00
      NBGGL  = 0                                                          40.51
      FSHAPE = 2                                                          40.00
      DSHAPE = 2                                                          40.00
      PSHAPE(1) = 3.3                                                     40.00
      PSHAPE(2) = 0.1                                                     40.00
!     ***** input grids *****
      DO 80 IGRID = 1, NUMGRD
        XPG(IGRID)    = 0.
        YPG(IGRID)    = 0.
        ALPG(IGRID)   = 0.
        COSPG(IGRID)  = 1.                                                40.13
        SINPG(IGRID)  = 0.
        DXG(IGRID)    = 0.
        DYG(IGRID)    = 0.
        MXG(IGRID)    = 0
        MYG(IGRID)    = 0
        LEDS(IGRID)   = 0
        STAGX(IGRID)  = 0.                                                30.21
        STAGY(IGRID)  = 0.                                                30.21
        EXCFLD(IGRID) = -1.E20                                            30.60
        IFLDYN(IGRID) = 0                                                 40.00
        IFLTIM(IGRID) = -1.E20                                            40.00
  80  CONTINUE
!     ***** computational grid *****
      OPTG   = 1                                                          40.80
      MXC    = 0
      MYC    = 0
      MXCGL  = 0                                                          40.30
      MYCGL  = 0                                                          40.30
      MXF    = 1                                                          40.30
      MXL    = 0                                                          40.80 40.30
      MYF    = 1                                                          40.30
      MYL    = 0                                                          40.80 40.30
      MSC    = 0
      MDC    = 0
      MTC    = 1
      ICOMP  = 1
      ALPC   = 0.
      FULCIR = .TRUE.
      SPDIR1 = 0.
!     number of points needed in computational stencil:
      ICMAX  = 3
!     ***** numerical scheme *****
      NCOR   = 1
      NSTATM = -1                                                         40.00
      NSTATC = -1                                                         40.00
      NCOMPT = 0                                                          40.41
!
!     initialise number of iterations stationary and nonstationary        40.03
      MXITST = 50                                                         40.80 40.03
      MXITNS = 1                                                          40.03
      ITERMX = MXITST                                                     40.03
      ICUR   = 0
      IDIF   = 0
      IINC   = 0
!
!     --- meaning IREFR:                                                  40.31
!         IREFR = -1: limiter on Ctheta activated                         30.80
!         IREFR =  1: No limiter on Ctheta                                30.80
!         IREFR =  0: No refraction                                       30.80
!
      IREFR  = 1                                                          40.02
      ITFRE  = 1
      IWIND  = 0
      IGEN   = 3                                                          32.06
      IQUAD  = 2                                                          30.72
      IWCAP  = 1                                                          21/MAY
      ISURF  = 1
      IBOT   = 0
      ITRIAD = 0
      IVEG   = 0                                                          40.55
      VARWI  = .FALSE.
      VARFR  = .FALSE.
      VARWLV = .FALSE.                                                    20.38
      VARAST = .FALSE.       ! True means spatially variable air-sea t.d. 40.03
      VARNPL = .FALSE.                                                    40.55
      U10    = 0.
      WDIP   = 0.
      INRHOG = 0                                                          30.20
      DEPMIN = 0.05
      SY0    = 3.3
      SIGMAG = 0.1
      XOFFS  = 0.
      YOFFS  = 0.
      LXOFFS = .FALSE.
      DYNDEP = .FALSE.                                                    40.00
      LWDATE = 0                                                          41.13
      MXOUTAR = 0
!
      FBS%NBS = -999                                                      40.31
!
!     Set the defaults for the MDIA:                                      40.17

      MDIA  = 6                                                           40.17
      IF(ALLOCATED(LAMBDA)) DEALLOCATE(LAMBDA)
      IF(ALLOCATED(CNL4_1)) DEALLOCATE(CNL4_1)
      IF(ALLOCATED(CNL4_2)) DEALLOCATE(CNL4_2)
      ALLOCATE(LAMBDA(MDIA),CNL4_1(MDIA),CNL4_2(MDIA))                    40.17
      LAMBDA = (/0.08,0.09,0.11,0.15,0.16,0.29/)                          40.17
      CNL4_1 = (/8.77,-13.82,10.02,-15.92,14.41,0.65/)                    40.17
      CNL4_2 = CNL4_1                                                     40.17
      CNL4_1 = CNL4_1 * ((2.*PI)**9)                                      40.17
      CNL4_2 = CNL4_2 * ((2.*PI)**9)                                      40.17
!
!     *** Initial conditions ***
      ICOND = 0                                                          060697
!
      BNAUT  = .FALSE.                                                    32.01
      BNDCHK = .TRUE.                                                     32.01
      BRESCL = .TRUE.                                                     40.00
      ONED   = .FALSE.                                                    32.02
      ACUPDA = .TRUE.                                                     40.07
      OFFSRC = .FALSE.                                                    40.80
      LADDS  = .FALSE.                                                    40.85
      HSRERR = 0.1                                                        32.01
!
!     higher order propagation and spherical coordinates                  33.08
!
      PROJ_METHOD = 0                                                     33.09
      PROPSS = 2                                                          33.08
      PROPSN = 3                                                          33.08
      PROPSC = 1                                                          33.08
      PROPSL = 1                                                          33.08
      PROPFL = 0                                                          40.23
      WAVAGE = 0.                                                         33.08
      KSPHER = 0                                                          33.09
      KREPTX = 0                                                          33.09
      REARTH = 2.E7/PI                                                    33.09
      LENDEG = 2.E7/180.                                                  33.09
!
!     *** setup flag ***                                                  32.02
      LSETUP = 0                                                          32.02
!
!     *** flag for setup convergence                                      30.82
!
      CSETUP = .TRUE.                                                     30.82
!
!     flag for frequency dependent surf breaking                          41.06
!
      IFRSRF = 0                                                          41.06
!
!     PSETUP(1) is currently unused, but can be used as setup nesting flag
!     PSETUP(2) is the user defined correction for the level of the setup
!
      PSETUP(1) = 0.0                                                     30.82
      PSETUP(2) = 0.0                                                     30.82
!
!     *** ACCURACY criterion ***
!
!     *** relative error in significant wave height and mean period ***
      PNUMS(1)  = 0.02                                                    30.82
!     *** absolute error in significant wave heigth (m) ***
      PNUMS(2)  = 0.03
!     *** absolute error in mean wave period (s) ***
      PNUMS(3)  = 0.3
!     *** total number of wet gridpoints were accuracy has ***
!     *** been reached                                     ***
      PNUMS(4)  = 98.00
!
!     *** DIFFUSION schemes ***
!
!     *** Numerical diffusion over theta ***
      PNUMS(6)  = 0.5                                                     20.78
!     *** Numerical diffusion over sigma ***
      PNUMS(7)  = 0.5                                                     20.78
!     *** Explicit or implicit scheme in frequency space ***
!     *** default = implicit : PNUMS(8) = 1              ***
      PNUMS(8) = 1.
!     *** diffusion coefficient for explicit scheme ***
      PNUMS(9) = 0.01
!
!     *** parameters for the SIP solver                        ***
!
!     *** Required accuracy to terminate the solver            ***
!     ***                                                      ***
!     ***  || Ax-b ||  <  eps2 * || b ||                       ***
!     ***                                                      ***
!     ***  eps2 = PNUMS(12)                                    ***
!
!     *** PNUMS(13) output for the solver. Possible values:    ***
!     ***     <0  : no output                                  ***
!     ***      0  : only fatal errors will be printed          ***
!     ***      1  : additional information about the iteration ***
!     ***           is printed                                 ***
!     ***      2  : gives a maximal amount of output           ***
!     ***           concerning the iteration process           ***
!
!     *** PNUMS(14) : maximum number of iterations             ***
!
      PNUMS(12) = 1.E-4
      PNUMS(13) = 0.
      PNUMS(14) = 20.
!
!     For the setup calculation, next parameters for the solver are used:
!
!     PNUMS(23) : required accuracy to terminate the solver
!     PNUMS(24) : output for the solver (see PNUMS(13) for meanings)
!     PNUMS(25) : maximum number of iterations
!
      PNUMS(23) = 1.E-6                                                   40.41 30.82
      PNUMS(24) = 0.                                                      30.82
      PNUMS(25) = 1000.                                                   40.41 30.82
!
!     Maximum growth in spectral bin
!     The value is the default in the command GEN3 KOM
!
      PNUMS(20) = 0.1                                                     30.72
!
!     Added coefficient for use with limiter on action (Qb switch)        40.16
!
      PNUMS(28) = 1.                                                      40.23
!
!     *** set the values of PNUMS that are not used equal 0. ***
!
      PNUMS(5)  = 0.
!
!     The allowed global errors in the iteration procedure:               30.82
!     PNUMS(15) for Hs and PNUMS(16) for Tm01                             30.82
!
      PNUMS(15) = 0.02                                                    30.82
      PNUMS(16) = 0.02                                                    30.82
!
!     coefficient for limitation of Ctheta                                30.80
!     default no limitation on refraction                                 40.02
!
      PNUMS(17) = -1.                                                     40.02
!
!     Limitation on Froude number; current velocity is reduced if greater
!     than Pnums(18)*Sqrt(grav*depth)
!
      PNUMS(18) = 0.8                                                     30.50
!
!     *** CFL criterion for explicit scheme in frequency space ***
!
      PNUMS(19) = 0.5 * sqrt (2.)
!
!     --- coefficient for type stopping criterion                         40.41
!
      PNUMS(21) = -9.                                                     40.80 40.41
!
!     --- under-relaxation factor
!
      PNUMS(30) = 0.00                                                    40.23
!
!     --- parameters for limiting Ctheta
!
      PNUMS(26) = 0.2                                                     41.06
      PNUMS(27) = 2.0                                                     41.06
      PNUMS(29) = 0.0                                                     41.06
!
!     --- computation of Ctheta based on wave number
!
      PNUMS(32) = 0.                                                      41.07
!
!     *** (1) and (2): Komen et al. (1984) formulation ***
!
      PWCAP(1)  = 2.36E-5
      PWCAP(2)  = 3.02E-3
      PWCAP(9)  = 2.                                                      34.00
      PWCAP(10) = 0.                                                      34.00
      PWCAP(11) = 1.                                                      34.00
!
!     *** (3): Coefficient for Janssen(1989,1991) formulation ***
!     ** according to Komen et al. (1994) ***
!
      PWCAP(3)  = 4.5
      PWCAP(4)  = 0.5
!
!     *** (5): Coefficient for Longuet-Higgins ***
!
      PWCAP(5) = 1.
!
!     *** (6): ALPHA in Battjes/Janssen ***
!
      PWCAP(6) = 0.88
      PWCAP(7) = 1.
      PWCAP(8) = 0.75
!
!     PBOT(1)   = 0.005          modified
      PBOT(1)   = 0.0                                                     20.68
      PBOT(2)   = 0.015                                                   20.68
      PBOT(3)   = 0.067
      PBOT(4)   = -0.08
      PBOT(5)   = 0.05
!
      PSURF(1)  = 1.0
      PSURF(2)  = 0.73                                                    20.67
!
!     triad interactions
      PTRIAD(1)  = 0.05                                                   40.61 30.82
      PTRIAD(2)  = 2.5                                                    40.61 40.56 30.82
      PTRIAD(3)  = 10.                                                    40.23
      PTRIAD(4)  = 0.2                                                    40.13
      PTRIAD(5)  = 0.01                                                   40.23
!
!     quadruplet interactions
      PQUAD(1) = 0.25                                                     34.00
      PQUAD(2) = 3.E7                                                     34.00
      PQUAD(3) = 5.5                                                      34.00
      PQUAD(4) = 0.833                                                    34.00
      PQUAD(5) = -1.25                                                    34.00
!
      PWIND(1)  = 188.0
      PWIND(2)  = 0.59
      PWIND(3)  = 0.12
      PWIND(4)  = 250.0
      PWIND(5)  = 0.0023
      PWIND(6)  = -0.223
      PWIND(7)  = 0.
      PWIND(8)  = -0.56
      PWIND(10) = 0.0036
      PWIND(11) = 0.00123
      PWIND(12) = 1.0
      PWIND(13) = 0.13
!     *** Janssen (1991) wave growth model ***
!     *** alpha ***
      PWIND(14) = 0.01
!      PWIND(14) = 0.0144
!     *** Charnock: Von Karman constant ***
      PWIND(15) = 0.41
!     *** rho air (density) ****
      PWIND(16) = 1.28
!     *** rho water (density) ***
      PWIND(17) = RHO
      PWIND(9)  = PWIND(16) / RHO
!     Coefficient in front of A term in 3d gen. growth term
!     default is 0; can be made non-zero in command GEN3 or GROWTH
      PWIND(31) = 0.                                                      7/MAR
!
!     --- coefficients for diffraction approximation                      40.21
!
      IDIFFR    = 0                                                       40.21
      PDIFFR(:) = 0.                                                      40.21
!
!     pointers in array COMPDA
      JDISS = 2
      JUBOT = 3
      JQB   = 4
      JSTP  = 5
      JDHS  = 6
      JDP1  = 7
      JDP2  = 8
      JVX1  = 9
      JVY1  = 10
      JVX2  = 11
      JVY2  = 12
      JVX3  = 13
      JVY3  = 14
      JDP3  = 15
      JWX2  = 16                                                          40.00
      JWY2  = 17                                                          40.00
      JWX3  = 18                                                          40.00
      JWY3  = 19                                                          40.00
      JDTM  = 20                                                          40.00
      JLEAK = 21                                                          40.00
      JWLV1 = 22                                                          40.00
      JWLV3 = 23                                                          40.00
      JWLV2 = 24                                                          40.00
      JHSIBC = 25                                                         40.00
      JHS    = 26                                                         40.13
      JURSEL = 27                                                         40.13
      MCMVAR = 27                                                         40.65 40.61 40.51 40.13
!     subarray sequence number 1 is used only for unused subarrays        40.13
      JFRC2 = 1                                                           40.00
      JFRC3 = 1                                                           40.00
      JUSTAR= 1                                                           30.22
      JZEL  = 1                                                           30.22
      JTAUW = 1                                                           30.22
      JCDRAG= 1                                                           30.22
!     added for air-sea temp. diff.:
      JASTD2= 1                                                           40.03
      JASTD3= 1                                                           40.03
!     added for number of plants per square meter                         40.55
      JNPLA2= 1                                                           40.55
      JNPLA3= 1                                                           40.55
!
!     *** added for wave setup ***                                        32.01
!
      JSETUP = 1                                                          32.02
      JDPSAV = 1                                                          32.02
!
!     --- added for output purposes                                       40.65
!
      JPBOT  = 1                                                          40.65 40.51
      JBOTLV = 1                                                          40.65
      JDSXB  = 1                                                          40.65 40.61
      JDSXS  = 1                                                          40.65 40.61
      JDSXW  = 1                                                          40.65 40.61
      JDSXV  = 1                                                          40.65 40.61
      JGENR  = 1                                                          40.85
      JGSXW  = 1                                                          40.85
      JREDS  = 1                                                          40.85
      JRSXQ  = 1                                                          40.85
      JRSXT  = 1                                                          40.85
      JTRAN  = 1                                                          40.85
      JTSXG  = 1                                                          40.85
      JTSXT  = 1                                                          40.85
      JTSXS  = 1                                                          40.85
      JRADS  = 1                                                          40.85
!
!     Next pointers are for the plot of the source terms (SWTSDA)         40.31
!
      JPWNDA = 1
      JPWNDB = 2
      JPWCAP = 3
      JPBTFR = 4
      JPWBRK = 5
      JP4S   = 6
      JP4D   = 7
      JPTRI  = 8
      JPVEGT = 9                                                          40.55
      MTSVAR = 9
!
!     ***** test output control *****
      ITEST  = 1                                                          40.41
      INTES  = 0                                                          30.50
      ICOTES = 0                                                          30.50
      IOUTES = 0                                                          30.50
      LTRACE = .FALSE.
      TESTFL = .FALSE.
      NPTST  = 0
      NPTSTA = 1
      LXDMP  = -1
      LYDMP  = 0
      NEGMES = 0
      MAXMES = 200
      IFPAR = 0                                                           40.00
      IFS1D = 0                                                           40.00
      IFS2D = 0                                                           40.00
!     number of obstacles initialised at 0                                40.13
      NUMOBS = 0                                                          40.13
!     ***** output *****
      IUBOTR = 0
!     ***** plot output *****
!
      DO IVT = 1, NMOVAR
        OVKEYW(IVT) = 'XXXX'                                              40.00
      ENDDO
!
!     properties of output variables
!
      IVTYPE = 1
!     keyword used in SWAN command
      OVKEYW(IVTYPE) = 'XP'                                               40.00
!     short name
      OVSNAM(IVTYPE) = 'Xp'
!     long name
      OVLNAM(IVTYPE) = 'X user coordinate'
!     unit name
      OVUNIT(IVTYPE) = UL
!     type (scalar/vector etc.)
      OVSVTY(IVTYPE) = 1
!     lower and upper limit
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
!     lowest and highest expected value
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
!     exception value
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 2
      OVKEYW(IVTYPE) = 'YP'                                               40.00
      OVSNAM(IVTYPE) = 'Yp'
      OVLNAM(IVTYPE) = 'Y user coordinate'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 3
      OVKEYW(IVTYPE) = 'DIST'                                             40.00
      OVSNAM(IVTYPE) = 'Dist'
      OVLNAM(IVTYPE) = 'distance along output curve'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 4
      OVKEYW(IVTYPE) = 'DEP'                                              40.00
      OVSNAM(IVTYPE) = 'Depth'
      OVLNAM(IVTYPE) = 'Depth'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 5
      OVKEYW(IVTYPE) = 'VEL'                                              40.00
      OVSNAM(IVTYPE) = 'Vel'
      OVLNAM(IVTYPE) = 'Current velocity'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -2.
      OVHEXP(IVTYPE) = 2.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 6
      OVKEYW(IVTYPE) = 'UBOT'                                             40.00
      OVSNAM(IVTYPE) = 'Ubot'
      OVLNAM(IVTYPE) = 'RMS of orbital velocity amplitude at the bottom'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -10.
!
      IVTYPE = 7
      OVKEYW(IVTYPE) = 'DISS'                                             40.00
      OVSNAM(IVTYPE) = 'Dissip'
      OVLNAM(IVTYPE) = 'Energy dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 8
      OVKEYW(IVTYPE) = 'QB'                                               40.00
      OVSNAM(IVTYPE) = 'Qb'
      OVLNAM(IVTYPE) = 'Fraction breaking waves'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -1.
!
      IVTYPE = 9
      OVKEYW(IVTYPE) = 'LEA'                                              40.00
      OVSNAM(IVTYPE) = 'Leak'
      OVLNAM(IVTYPE) = 'Energy leak over spectral boundaries'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 10
      OVKEYW(IVTYPE) = 'HS'                                               40.00
      OVSNAM(IVTYPE) = 'Hsig'
      OVLNAM(IVTYPE) = 'Significant wave height'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = -9.
!                                                                         modified 10.09
      IVTYPE = 11
      OVKEYW(IVTYPE) = 'TM01'                                             40.00
      OVSNAM(IVTYPE) = 'Tm01'                                             20.81
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 12
      OVKEYW(IVTYPE) = 'RTP'                                              40.00
      OVSNAM(IVTYPE) = 'RTpeak'                                           40.41 20.81
      OVLNAM(IVTYPE) = 'Relative peak period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 13
      OVKEYW(IVTYPE) = 'DIR'                                              40.00
      OVSNAM(IVTYPE) = 'Dir'
      OVLNAM(IVTYPE) = 'Average wave direction'
      OVUNIT(IVTYPE) = UDI                                                30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 14
      OVKEYW(IVTYPE) = 'PDI'                                              40.00
      OVSNAM(IVTYPE) = 'PkDir'
      OVLNAM(IVTYPE) = 'direction of the peak of the spectrum'
      OVUNIT(IVTYPE) = UDI                                                30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 15
      OVKEYW(IVTYPE) = 'TDI'                                              40.00
      OVSNAM(IVTYPE) = 'TDir'
      OVLNAM(IVTYPE) = 'direction of the energy transport'
      OVUNIT(IVTYPE) = UDI                                                30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 16
      OVKEYW(IVTYPE) = 'DSPR'                                             40.00
      OVSNAM(IVTYPE) = 'Dspr'
      OVLNAM(IVTYPE) = 'directional spreading'
      OVUNIT(IVTYPE) = UDI                                                30.72
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 60.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 17
      OVKEYW(IVTYPE) = 'WLEN'                                             40.00
      OVSNAM(IVTYPE) = 'Wlen'
      OVLNAM(IVTYPE) = 'Average wave length'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 200.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 18
      OVKEYW(IVTYPE) = 'STEE'                                             40.00
      OVSNAM(IVTYPE) = 'Steepn'
      OVLNAM(IVTYPE) = 'Wave steepness'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 19
      OVKEYW(IVTYPE) = 'TRA'                                              40.00
      OVSNAM(IVTYPE) = 'Transp'
      OVLNAM(IVTYPE) = 'Wave energy transport'
      OVUNIT(IVTYPE) = 'm3/s'                                             40.00
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 20
      OVKEYW(IVTYPE) = 'FOR'                                              40.00
      OVSNAM(IVTYPE) = 'WForce'
      OVLNAM(IVTYPE) = 'Wave driven force per unit surface'
      OVUNIT(IVTYPE) = UF                                                 30.72
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -1.E5
      OVULIM(IVTYPE) =  1.E5
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -9.                                                40.80
!
      IVTYPE = 21
      OVKEYW(IVTYPE) = 'AAAA'                                             40.00
      OVSNAM(IVTYPE) = 'AcDens'
      OVLNAM(IVTYPE) = 'spectral action density'
      OVUNIT(IVTYPE) = 'm2s'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 22
      OVKEYW(IVTYPE) = 'EEEE'                                             40.00
      OVSNAM(IVTYPE) = 'EnDens'
      OVLNAM(IVTYPE) = 'spectral energy density'
      OVUNIT(IVTYPE) = 'm2'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 23
      OVKEYW(IVTYPE) = 'AAAA'                                             40.00
      OVSNAM(IVTYPE) = 'Aux'
      OVLNAM(IVTYPE) = 'auxiliary variable'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 24
      OVKEYW(IVTYPE) = 'XC'                                               40.00
      OVSNAM(IVTYPE) = 'Xc'
      OVLNAM(IVTYPE) = 'X computational grid coordinate'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 25
      OVKEYW(IVTYPE) = 'YC'                                               40.00
      OVSNAM(IVTYPE) = 'Yc'
      OVLNAM(IVTYPE) = 'Y computational grid coordinate'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 26
      OVKEYW(IVTYPE) = 'WIND'                                             40.00
      OVSNAM(IVTYPE) = 'Windv'
      OVLNAM(IVTYPE) = 'Wind velocity at 10 m above sea level'
      OVUNIT(IVTYPE) = UV                                                 30.72
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -50.
      OVHEXP(IVTYPE) = 50.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 27
      OVKEYW(IVTYPE) = 'FRC'                                              40.00
      OVSNAM(IVTYPE) = 'FrCoef'
      OVLNAM(IVTYPE) = 'Bottom friction coefficient'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!                                                                     new 10.09
      IVTYPE = 28
      OVKEYW(IVTYPE) = 'RTM01'                                            40.00
      OVSNAM(IVTYPE) = 'RTm01'                                            20.81
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 29                                                         20.28
      OVKEYW(IVTYPE) = 'EEEE'                                             40.00
      OVSNAM(IVTYPE) = 'EnDens'
      OVLNAM(IVTYPE) = 'energy density integrated over direction'         20.28
      OVUNIT(IVTYPE) = 'm2'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 30                                                         20.52
      OVKEYW(IVTYPE) = 'DHS'                                              40.00
      OVSNAM(IVTYPE) = 'dHs'
      OVLNAM(IVTYPE) = 'difference in Hs between iterations'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 31                                                         20.52
      OVKEYW(IVTYPE) = 'DRTM01'                                           40.00
      OVSNAM(IVTYPE) = 'dTm'
      OVLNAM(IVTYPE) = 'difference in Tm between iterations'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 2.
      OVEXCV(IVTYPE) = -9.
!                                                                         20.61
      IVTYPE = 32
      OVKEYW(IVTYPE) = 'TM02'                                             40.00
      OVSNAM(IVTYPE) = 'Tm02'
      OVLNAM(IVTYPE) = 'Zero-crossing period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!                                                                         20.61
      IVTYPE = 33
      OVKEYW(IVTYPE) = 'FSPR'                                             40.00
      OVSNAM(IVTYPE) = 'FSpr'                                             20.67
      OVLNAM(IVTYPE) = 'Frequency spectral width (Kappa)'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 34                                                         20.67
      OVKEYW(IVTYPE) = 'URMS'                                             40.00
      OVSNAM(IVTYPE) = 'Urms'
      OVLNAM(IVTYPE) = 'RMS of orbital velocity at the bottom'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 35                                                         30.22
      OVKEYW(IVTYPE) = 'UFRI'                                             40.00
      OVSNAM(IVTYPE) = 'Ufric'
      OVLNAM(IVTYPE) = 'Friction velocity'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 36                                                         30.22
      OVKEYW(IVTYPE) = 'ZLEN'                                             40.00
      OVSNAM(IVTYPE) = 'Zlen'
      OVLNAM(IVTYPE) = 'Zero velocity thickness of boundary layer'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 37                                                         30.22
      OVKEYW(IVTYPE) = 'TAUW'                                             40.00
      OVSNAM(IVTYPE) = 'TauW'
      OVLNAM(IVTYPE) = '    '
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 38                                                         30.22
      OVKEYW(IVTYPE) = 'CDRAG'                                            40.00
      OVSNAM(IVTYPE) = 'Cdrag'
      OVLNAM(IVTYPE) = 'Drag coefficient'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
!     *** wave-induced setup ***                                          32.02
!
      IVTYPE = 39                                                         32.02
      OVKEYW(IVTYPE) = 'SETUP'                                            40.00
      OVSNAM(IVTYPE) = 'Setup'                                            32.02
      OVLNAM(IVTYPE) = 'Setup due to waves'                               32.02
      OVUNIT(IVTYPE) = 'm'                                                32.02
      OVSVTY(IVTYPE) = 1                                                  32.02
      OVLLIM(IVTYPE) = -1.                                                32.02
      OVULIM(IVTYPE) = 1.                                                 32.02
      OVLEXP(IVTYPE) = -1.                                                32.02
      OVHEXP(IVTYPE) = 1.                                                 32.02
      OVEXCV(IVTYPE) = -9.                                                32.02
!
      IVTYPE = 40                                                         40.00
      OVKEYW(IVTYPE) = 'TIME'                                             40.00
      OVSNAM(IVTYPE) = 'Time'
      OVLNAM(IVTYPE) = 'Date-time'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -99999.
!
      IVTYPE = 41                                                         40.00
      OVKEYW(IVTYPE) = 'TSEC'                                             40.00
      OVSNAM(IVTYPE) = 'Tsec'
      OVLNAM(IVTYPE) = 'Time in seconds from reference time'
      OVUNIT(IVTYPE) = 's'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100000.
      OVLEXP(IVTYPE) = -100000.
      OVHEXP(IVTYPE) = 1000000.
      OVEXCV(IVTYPE) = -99999.
!                                                        new              40.00
      IVTYPE = 42
      OVKEYW(IVTYPE) = 'PER'                                              40.00
      OVSNAM(IVTYPE) = 'Period'                                           40.00
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!                                                        new              40.00
      IVTYPE = 43
      OVKEYW(IVTYPE) = 'RPER'                                             40.00
      OVSNAM(IVTYPE) = 'RPeriod'                                          40.41 40.00
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 44                                                         40.00
      OVKEYW(IVTYPE) = 'HSWE'                                             40.00
      OVSNAM(IVTYPE) = 'Hswell'
      OVLNAM(IVTYPE) = 'Wave height of swell part'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 45
      OVKEYW(IVTYPE) = 'URSELL'                                           40.03
      OVSNAM(IVTYPE) = 'Ursell'
      OVLNAM(IVTYPE) = 'Ursell number'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 46
      OVKEYW(IVTYPE) = 'ASTD'                                             40.03
      OVSNAM(IVTYPE) = 'ASTD'
      OVLNAM(IVTYPE) = 'Air-Sea temperature difference'
      OVUNIT(IVTYPE) = 'K'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -50.
      OVULIM(IVTYPE) =  50.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 47                                                         40.41
      OVKEYW(IVTYPE) = 'TMM10'
      OVSNAM(IVTYPE) = 'Tm_10'
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 48                                                         40.41
      OVKEYW(IVTYPE) = 'RTMM10'
      OVSNAM(IVTYPE) = 'RTm_10'
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 49
      OVKEYW(IVTYPE) = 'DIFPAR'                                           40.21
      OVSNAM(IVTYPE) = 'DifPar'
      OVLNAM(IVTYPE) = 'Diffraction parameter'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -50.
      OVULIM(IVTYPE) =  50.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 50                                                         40.51
      OVKEYW(IVTYPE) = 'TMBOT'
      OVSNAM(IVTYPE) = 'TmBot'
      OVLNAM(IVTYPE) = 'Near bottom wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 51                                                         40.51
      OVKEYW(IVTYPE) = 'WATL'
      OVSNAM(IVTYPE) = 'Watlev'
      OVLNAM(IVTYPE) = 'Water level'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 52                                                         40.51
      OVKEYW(IVTYPE) = 'BOTL'
      OVSNAM(IVTYPE) = 'Botlev'
      OVLNAM(IVTYPE) = 'Bottom level'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 53                                                         40.51
      OVKEYW(IVTYPE) = 'TPS'
      OVSNAM(IVTYPE) = 'TPsmoo'
      OVLNAM(IVTYPE) = 'Relative peak period (smooth)'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 54
      OVKEYW(IVTYPE) = 'DISB'                                             40.61
      OVSNAM(IVTYPE) = 'Sfric'                                            40.85
      OVLNAM(IVTYPE) = 'Bottom friction dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 55
      OVKEYW(IVTYPE) = 'DISSU'                                            40.61
      OVSNAM(IVTYPE) = 'Ssurf'                                            40.85
      OVLNAM(IVTYPE) = 'Surf breaking dissipation'                        40.85
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 56
      OVKEYW(IVTYPE) = 'DISW'                                             40.61
      OVSNAM(IVTYPE) = 'Swcap'                                            40.85
      OVLNAM(IVTYPE) = 'Whitecapping dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 57
      OVKEYW(IVTYPE) = 'DISV'                                             40.61
      OVSNAM(IVTYPE) = 'Sveg'                                             40.85
      OVLNAM(IVTYPE) = 'Vegetation dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 58                                                         40.64
      OVKEYW(IVTYPE) = 'QP'
      OVSNAM(IVTYPE) = 'Qp'
      OVLNAM(IVTYPE) = 'Peakedness'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 59                                                         40.64
      OVKEYW(IVTYPE) = 'BFI'
      OVSNAM(IVTYPE) = 'BFI'
      OVLNAM(IVTYPE) = 'Benjamin-Feir index'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1000.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 60
      OVKEYW(IVTYPE) = 'GENE'                                             40.85
      OVSNAM(IVTYPE) = 'Genera'
      OVLNAM(IVTYPE) = 'Energy generation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 61
      OVKEYW(IVTYPE) = 'GENW'                                             40.85
      OVSNAM(IVTYPE) = 'Swind'
      OVLNAM(IVTYPE) = 'Wind source term'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 62
      OVKEYW(IVTYPE) = 'REDI'                                             40.85
      OVSNAM(IVTYPE) = 'Redist'
      OVLNAM(IVTYPE) = 'Energy redistribution'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 63
      OVKEYW(IVTYPE) = 'REDQ'                                             40.85
      OVSNAM(IVTYPE) = 'Snl4'
      OVLNAM(IVTYPE) = 'Total absolute 4-wave interaction'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 64
      OVKEYW(IVTYPE) = 'REDT'                                             40.85
      OVSNAM(IVTYPE) = 'Snl3'
      OVLNAM(IVTYPE) = 'Total absolute 3-wave interaction'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 65
      OVKEYW(IVTYPE) = 'PROPA'                                            40.85
      OVSNAM(IVTYPE) = 'Propag'
      OVLNAM(IVTYPE) = 'Energy propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 66
      OVKEYW(IVTYPE) = 'PROPX'                                            40.85
      OVSNAM(IVTYPE) = 'Propxy'
      OVLNAM(IVTYPE) = 'xy-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 67
      OVKEYW(IVTYPE) = 'PROPT'                                            40.85
      OVSNAM(IVTYPE) = 'Propth'
      OVLNAM(IVTYPE) = 'theta-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 68
      OVKEYW(IVTYPE) = 'PROPS'                                            40.85
      OVSNAM(IVTYPE) = 'Propsi'
      OVLNAM(IVTYPE) = 'sigma-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 69
      OVKEYW(IVTYPE) = 'RADS'                                             40.85
      OVSNAM(IVTYPE) = 'Radstr'
      OVLNAM(IVTYPE) = 'Radiation stress'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 70
      OVKEYW(IVTYPE) = 'NPL'                                            ! 41.12
      OVSNAM(IVTYPE) = 'Nplant'
      OVLNAM(IVTYPE) = 'Plants per m2'
      OVUNIT(IVTYPE) = '1/m2'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 71
      OVKEYW(IVTYPE) = 'LWAVP'                                            41.15
      OVSNAM(IVTYPE) = 'Lwavp'
      OVLNAM(IVTYPE) = 'Peak wave length'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 200.
      OVEXCV(IVTYPE) = -9.
!
!     various parameters for computation of output quantities             40.00
!
!     reference time for TSEC
      OUTPAR(1) = 0.
!     power in expression for PER and RPER
!     previous name: SPCPOW
      OUTPAR(2) = 1.
!     power in expression for WLEN
!     previous name: AKPOWR
      OUTPAR(3) = 1.
!     indicator for direction
!     =0: direction always w.r.t. user coordinates; =1: dir w.r.t. frame
      OUTPAR(4) = 0.
!     frequency limit for swell
      OUTPAR(5) = 0.1
!     0=integration over [0,inf], 1=integration over [fmin,fmax]
      OUTPAR(6:20) = 0.                                                   40.87
!     lower bound of integration range for output parameters
      OUTPAR(21:35) = 0.                                                  40.87
!     upper bound of integration range for output parameters
      OUTPAR(36:50) = 1000.                                               40.87
!
      RETURN
! * end of subroutine SWINIT *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID ,YCGRID ,        40.31
     &                    KGRPNT, KGRBND, SPCDIR, SPCSIG )                40.31
!                                                                      *
!***********************************************************************
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC                                                       40.31
      USE M_PARALL                                                        40.31
      USE M_DIFFR                                                         40.21
      USE SwanGriddata                                                    40.80

      INTEGER, INTENT(INOUT) :: KGRBND(*)                                 40.02
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
      REAL    SPCDIR(MDC,6)  ,    SPCSIG(MSC)                             40.31

      INTEGER   KGRPNT(MXC,MYC), CROSS(2,MCGRD)
      INTEGER   BGRIDP(6*NBGRPT)                                          40.31
      REAL      BSPECS(MDC,MSC,NBSPEC,2)                                  40.31
      INTEGER, ALLOCATABLE :: CROSSGL(:,:)                                40.31

      TYPE(BSDAT) , POINTER :: CURRBS                                     40.31
      TYPE(BGPDAT), POINTER :: CURBGP                                     40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT, 'SWPREP')

      IF ( ITEST.GE.100 ) THEN                                            40.31
         WRITE (PRTEST,*) ' BNAUT after SWREAD = ',BNAUT                  40.31
      END IF

!     coefficients for transformation from user coordinates to comp. coord.
      COSPC = COS(ALPC)
      SINPC = SIN(ALPC)
      XCP   = -XPC*COSPC - YPC*SINPC
      YCP   =  XPC*SINPC - YPC*COSPC
      ALCP  = -ALPC

!     wind direction w.r.t. computational grid
!
!     ALTMP = (WDIP - ALPC) / PI2              removed 30.50
      ALTMP = (WDIP) / PI2                                                13/JAN
      WDIC  = PI2 * (ALTMP - NINT(ALTMP))
!
      IF (MXC.LE.0 .AND. OPTG.NE.5) CALL MSGERR                           40.80
     & (3, 'no valid computational grid; check command CGRID')            40.80
      IF (MCGRD.LE.1 .AND. nverts.LE.0) CALL MSGERR                       40.80
     & (3, 'no valid comp. grid; check command READ BOT or READ UNSTRU')  40.80 30.50
!
      IF (LEDS(1).EQ.0) CALL MSGERR (3,'Bottom grid not defined')
      IF (LEDS(1).EQ.1) CALL MSGERR (3,'No bottom levels read')
      IF (IUBOTR.EQ.1 .AND. IBOT.EQ.0)
     &    CALL MSGERR (1,'Bottom friction not on, UBOT not computed')
!
      IF (LEDS(2).EQ.2) THEN
        IF (LEDS(3).NE.2)
     &  CALL MSGERR (3, 'VY not read, while VX is read')
!       ALBC  = ALPC - ALPG(2)
        ALBC  = - ALPG(2)
        COSVC = COS(ALBC)
        SINVC = SIN(ALBC)
      ENDIF
!
      IF (LEDS(4).EQ.2) VARFR = .TRUE.
!
      IF (LEDS(5).EQ.2) THEN
        IF (LEDS(6).NE.2)
     &  CALL MSGERR (3, 'WY not read, while WX is read')
        VARWI = .TRUE.
!       ALBC  = ALPC - ALPG(5)
        ALBC  = - ALPG(5)
        COSWC = COS(ALBC)
        SINWC = SIN(ALBC)
      ENDIF
!
      IF (LEDS(7).EQ.2) VARWLV = .TRUE.                                   20.38
!
      IF (IFLDYN(1).EQ.1 .OR. IFLDYN(7).EQ.1) DYNDEP = .TRUE.             40.00
!
      IF (OPTG.EQ.5) BNDCHK = .FALSE.                                     40.80
!
!     close input files containing stationary input fields                40.00
!
      DO IFLD = 1, NUMGRD
        IF (IFLDYN(IFLD).EQ.0 .AND. IFLNDS(IFLD).NE.0) THEN               40.00
          CLOSE (IFLNDS(IFLD))
          IFLNDS(IFLD) = 0
        ENDIF
      ENDDO
!
!     computation of tail factors for moments of action spectrum
!     IP=0: action int. IP=1: energy int. IP=2: first moment of energy etc.
!
      DO IP = 0, 3
        PPTAIL = PWTAIL(1) - REAL(IP)
        PWTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.71
      ENDDO
!
!     check location of output areas
!
      CALL SPRCON ( XCGRID, YCGRID, KGRPNT, KGRBND )                      40.31
!
!     *** find obstacles crossing the points in the stencil  ***
!     *** in structured grids                                ***          40.80
!
      IF (NUMOBS .GT. 0 .AND. OPTG.NE.5) THEN                             40.80
        DO INDX = 1, MCGRD                                               040697
          CROSS(1,INDX) = 0                                              040697
          CROSS(2,INDX) = 0                                              040697
        ENDDO                                                            040697
!
        ITMP1  = MXC                                                      40.31
        ITMP2  = MYC                                                      40.31
        ITMP3  = MCGRD                                                    40.31
        MXC    = MXCGL                                                    40.31
        MYC    = MYCGL                                                    40.31
        MCGRD  = MCGRDGL                                                  40.31
        ALLOCATE(CROSSGL(2,MCGRDGL))                                      40.31
        CROSSGL = 0                                                       40.31
        CALL OBSTMOVE (XGRDGL, YGRDGL, KGRPGL)                            40.31 40.09
        CALL SWOBST (XGRDGL, YGRDGL, KGRPGL, CROSSGL)                     40.31 30.7N
        MXC    = ITMP1                                                    40.31
        MYC    = ITMP2                                                    40.31
        MCGRD  = ITMP3                                                    40.31
!
        II = 1                                                            40.31
        DO IX = MXF, MXL                                                  40.31
           DO IY = MYF, MYL                                               40.31
              INDX = KGRPGL(IX,IY)                                        40.31
              IF ( INDX.NE.1 ) THEN                                       40.31
                 II = II + 1                                              40.31
                 CROSS(1:2,II) = CROSSGL(1:2,INDX)                        40.31
              END IF                                                      40.31
           END DO                                                         40.31
        END DO                                                            40.31
        DEALLOCATE(CROSSGL)                                               40.31
!
        IF (ITEST .GE. 120) THEN                                          060697
          WRITE(PRINTF,102)
 102      FORMAT('Links with obstacles crossing',/,
     &    '     COMP COORD           LINK         VALUE')
          DO IIYY =2,MYC
            DO IIXX = 2,MXC
            I1 = KGRPNT(IIXX,IIYY)
              DO I2 = 1,2
                IF (CROSS(I2,I1) .NE. 0)
     &          WRITE(PRINTF,101) IIXX,IIYY,I2,I1,CROSS(I2,I1)
 101            FORMAT(' POINT(',I4,',',I4,')',
     &          '  CROSS(',I3,',',I5,')  = ',I5)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
!
      CURRBS => FBS                                                       40.31
      DO                                                                  40.31
         NBS = CURRBS%NBS                                                 40.31
         IF (NBS.EQ.-999) EXIT                                            40.31
         SPPARM(1:4) = CURRBS%SPPARM(1:4)                                 40.31
         FSHAPE      = CURRBS%FSHAPE                                      40.31
         DSHAPE      = CURRBS%DSHAPE                                      40.31
         CALL SSHAPE (BSPECS(1,1,NBS,1), SPCSIG, SPCDIR,                  40.31
     &                FSHAPE, DSHAPE)                                     40.31
         IF (.NOT.ASSOCIATED(CURRBS%NEXTBS)) EXIT                         40.31
         CURRBS => CURRBS%NEXTBS                                          40.31
      END DO                                                              40.31

      BGRIDP = 0                                                          40.31
      IF (OPTG.NE.5) THEN                                                 40.80
         IF (NBGGL.EQ.0) NBGGL  = NBGRPT                                  40.51 40.31
         NBGRPT = 0                                                       40.31
         DO IX = MXF, MXL                                                 40.31
            DO IY = MYF, MYL                                              40.31
               INDX = KGRPGL(IX,IY)                                       40.31
               CURBGP => FBGP                                             40.31
               DO II = 1, NBGGL                                           40.31
                  INDXGR = CURBGP%BGP(1)                                  40.31
                  IF ( INDXGR.EQ.INDX ) THEN                              40.31
                     NBGRPT = NBGRPT + 1                                  40.31
                     BGRIDP(6*NBGRPT-5) = KGRPNT(IX-MXF+1,IY-MYF+1)       40.31
                     BGRIDP(6*NBGRPT-4) = CURBGP%BGP(2)                   40.31
                     BGRIDP(6*NBGRPT-3) = CURBGP%BGP(3)                   40.31
                     BGRIDP(6*NBGRPT-2) = CURBGP%BGP(4)                   40.31
                     BGRIDP(6*NBGRPT-1) = CURBGP%BGP(5)                   40.31
                     BGRIDP(6*NBGRPT  ) = CURBGP%BGP(6)                   40.31
                  END IF                                                  40.31
                  IF (.NOT.ASSOCIATED(CURBGP%NEXTBGP)) EXIT               40.31
                  CURBGP => CURBGP%NEXTBGP                                40.31
               END DO                                                     40.31
            END DO                                                        40.31
         END DO                                                           40.31
      ELSE                                                                40.80
         CURBGP => FBGP                                                   40.80
         DO II = 1, NBGRPT                                                40.80
            BGRIDP(6*II-5) = CURBGP%BGP(1)                                40.80
            BGRIDP(6*II-4) = CURBGP%BGP(2)                                40.80
            BGRIDP(6*II-3) = CURBGP%BGP(3)                                40.80
            BGRIDP(6*II-2) = CURBGP%BGP(4)                                40.80
            BGRIDP(6*II-1) = CURBGP%BGP(5)                                40.80
            BGRIDP(6*II  ) = CURBGP%BGP(6)                                40.80
            IF (.NOT.ASSOCIATED(CURBGP%NEXTBGP)) EXIT                     40.80
            CURBGP => CURBGP%NEXTBGP                                      40.80
         ENDDO                                                            40.80
      ENDIF                                                               40.80
!
!     --- allocate arrays for diffraction and set prop scheme to BSBT     40.41
!                                                                         40.21
      IF ( IDIFFR.EQ.1 ) THEN                                             40.21
         IF (.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(1:MCGRD))        40.21
         IF (.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(1:MCGRD))        40.21
         IF (.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(1:MCGRD))        40.21
         PROPSN = 1                                                       40.41
         PROPSS = 1                                                       40.41
      ELSE                                                                40.41
         IF (.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(0))              40.21
         IF (.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(0))              40.21
         IF (.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(0))              40.21
      END IF                                                              40.21
!
      RETURN
! * end of subroutine SWPREP *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SPRCON (XCGRID, YCGRID, KGRPNT, KGRBND)                  40.31 40.00
!                                                                      *
!***********************************************************************
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE OUTP_DATA                                                       40.31
      USE SwanGriddata                                                    40.80

      INTEGER :: KGRBND(*)                                                40.02
!
      REAL    :: XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      INTEGER I, J
      REAL    XR, YR
!
      INTEGER  KGRPNT(MXC,MYC)                                            40.02
      LOGICAL  SINBTG
      CHARACTER STYPE *1
      TYPE(OPSDAT), POINTER :: CUOPS                                      40.31
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT,'SPRCON')
!
!     ***** test location of computational grid *****
!
      IF (OPTG .EQ. 1) THEN                                               30.21
!
!       regular grid
!
        IF (ONED) THEN                                                    32.02
!
!         For 1D-version:
!         *** Check angles of bottom grid and computational grid ***      32.02
!
          IF ( ABS((ALPG(1) - ALPC)) .GT. 0.0017 ) THEN                   32.02
            CALL MSGERR( 2, ' Difference between angle of bottom grid'//  32.02
     &                      ' (alpinp) and computational grid (alpc)')    32.02
            CALL MSGERR( 2, ' greater than 0.1 degrees.')                 32.02
          ENDIF                                                           32.02
!
!         *** Check location of computational grid ***                    32.02
!
          DO  I=0,1                                                       30.72
            XR = I*XCLEN
            XP = XPC + XR*COSPC
            YP = YPC + XR*SINPC
            IF (.NOT.SINBTG(XP,YP) ) THEN
              CALL MSGERR(1,'Corner of comp grid outside bottom grid')
              WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS
            ENDIF                                                         30.72
          ENDDO
        ELSE                                                              32.02
!
!         two-dimensional case
!
          DO 11 I=0,1                                                     30.72
            DO 10 J=0,1
              XR = I*XCLEN
              YR = J*YCLEN
              XP = XPC + XR*COSPC - YR*SINPC
              YP = YPC + XR*SINPC + YR*COSPC
              IF (.NOT.SINBTG(XP,YP) ) THEN
                CALL MSGERR(1,'Corner of comp grid outside bottom grid')
                WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS
 6010           FORMAT (' Coordinates :',2F10.2)
              ENDIF                                                       30.72
   10       CONTINUE                                                      30.72
   11     CONTINUE
        ENDIF
!
        XR = 0.5*XCLEN
        YR = 0.5*YCLEN
        XP = XPC + XR*COSPC - YR*SINPC
        YP = YPC + XR*SINPC + YR*COSPC
        IF (.NOT. SINBTG(XP,YP) ) THEN
          CALL MSGERR (2,' Centre of comp. grid outside bottom grid')
        ENDIF
      ELSEIF (OPTG.EQ.5) THEN                                             40.80
!
!       --- check location of computational grid                          40.80
!
        DO I = 1, nverts                                                  40.80
           IF ( vmark(I) /= 0 ) THEN                                      40.80
              XP = xcugrd(I)                                              40.80
              YP = ycugrd(I)                                              40.80
              IF (.NOT.SINBTG(XP,YP) ) THEN                               40.80
                CALL MSGERR(1,'Corner of comp grid outside bottom grid')  40.80
                WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS                   40.80
              ENDIF                                                       40.80
           ENDIF                                                          40.80
        ENDDO                                                             40.80
      ENDIF                                                               30.21
!
!     ***** test location of output pointsets *****
      IF (LOPS) THEN                                                      40.31
        CUOPS => FOPS                                                     40.31
        DO                                                                40.31
!         --- get name of point set                                       40.31
          SNAME = CUOPS%PSNAME                                            40.31
          IF (ITEST.GE.80) WRITE (PRTEST, 12) SNAME                       40.31
  12      FORMAT (' test SPRCON ', A8)                                    40.31
!
!         bottom grid is excluded from the test
          IF (SNAME.EQ.'BOTTGRID') GOTO 100
!
!         computational grid is excluded from the test
          IF (SNAME.EQ.'COMPGRID') GOTO 100
!
!         wind grid is excluded from the test
          IF (SNAME.EQ.'WXGRID' .OR. SNAME.EQ.'WYGRID') GOTO 100
!
!         velocity grid is excluded from the test
          IF (SNAME.EQ.'VXGRID' .OR. SNAME.EQ.'VYGRID') GOTO 100
!
!         waterlevel grid is excluded from the test
          IF (SNAME.EQ.'WLEVGRID') GOTO 100
!
!         friction grid is excluded from the test
          IF (SNAME.EQ.'FRICGRID') GOTO 100
!
          STYPE = CUOPS%PSTYPE                                            40.31
!                                                                         32.02
!         *** Check other output locations ***                            32.02
!                                                                         32.02
          IF (STYPE.EQ.'F' .AND. OPTG .EQ. 1) THEN
!
!           check the four corners of the frame
!
            XQLEN = CUOPS%OPR(3)                                          40.31
            YQLEN = CUOPS%OPR(4)                                          40.31
            XPQ   = CUOPS%OPR(1)                                          40.31
            YPQ   = CUOPS%OPR(2)                                          40.31
            ALPQ  = CUOPS%OPR(5)                                          40.31
            COSPQ = COS(ALPQ)
            SINPQ = SIN(ALPQ)
            IF (ONED) THEN                                                32.02
              DO  I=0,1                                                   30.72
                XQ = I*XQLEN
                XP = XPQ + XQ*COSPQ
                YP = YPQ + XQ*SINPQ
                CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,       30.72
     &                         KGRBND)                                    40.00
              ENDDO                                                       30.72
            ELSE                                                          32.02
              DO 21 I=0,1                                                 30.72
                DO 20 J=0,1
                  XQ = I*XQLEN
                  YQ = J*YQLEN
                  XP = XPQ + XQ*COSPQ - YQ*SINPQ
                  YP = YPQ + XQ*SINPQ + YQ*COSPQ
                  CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,     30.72
     &                         KGRBND)                                    40.00
   20           CONTINUE                                                  30.72
   21         CONTINUE                                                    30.72
            ENDIF                                                         32.02
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'C') THEN
!
!           check first and last point of a curve
!
            MXK = CUOPS%MIP                                               40.31
            DO 30 IXK=1, MXK, MXK
              XP = CUOPS%XP(IXK)                                          40.31
              YP = CUOPS%YP(IXK)                                          40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,         30.72
     &                     KGRBND)                                        40.00
   30       CONTINUE
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'P' .OR. STYPE .EQ. 'U') THEN                    40.80
!
!           check all individual output points
!
            MIP = CUOPS%MIP                                               40.31
            DO 40 IP=1,MIP
              XP = CUOPS%XP(IP)                                           40.31
              YP = CUOPS%YP(IP)                                           40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,         30.72
     &                     KGRBND)                                        40.00
   40       CONTINUE
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'R') THEN
            MIP = CUOPS%MIP                                               40.31
            DO 50 IP=1,MIP, MIP
              XP = CUOPS%XP(IP)                                           40.31
              YP = CUOPS%YP(IP)                                           40.31
              XQ = CUOPS%XQ(IP)                                           40.31
              YQ = CUOPS%YQ(IP)                                           40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,         30.72
     &                     KGRBND)                                        40.00
              CALL SINUPT (SNAME, XQ, YQ, XCGRID, YCGRID, KGRPNT,         30.72
     &                     KGRBND)                                        40.00
   50       CONTINUE
          ENDIF
!         --------------------------------------------------------------
!         stype = 'N'     VERSION 20.63
!
          IF (STYPE.EQ.'N') THEN
            MIP   = CUOPS%MIP                                             40.31
            XQLEN = CUOPS%OPR(1)                                          40.31
            IF (XQLEN.NE.-999.) THEN                                      40.80
!              nested grid is regular, check corners
               YQLEN = CUOPS%OPR(2)                                       40.31
               XPQ   = CUOPS%OPR(3)                                       40.31
               YPQ   = CUOPS%OPR(4)                                       40.31
               ALPQ  = CUOPS%OPR(5)                                       40.31
               COSPQ = COS(ALPQ)
               SINPQ = SIN(ALPQ)
               DO I=0,1                                                   40.02
                 DO J=0,1                                                 40.02
                   XQ = I*XQLEN
                   YQ = J*YQLEN
                   XP = XPQ + XQ*COSPQ - YQ*SINPQ
                   YP = YPQ + XQ*SINPQ + YQ*COSPQ
                   CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,    30.72
     &                          KGRBND)                                   40.00
                 ENDDO                                                    40.02
               ENDDO                                                      40.02
            ELSE                                                          40.80
!              nested grid is unstructured, check whole outline
               DO IP=1,MIP                                                40.80
                  XP = CUOPS%XP(IP)                                       40.80
                  YP = CUOPS%YP(IP)                                       40.80
                  CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT,     40.80
     &                         KGRBND)                                    40.80
               ENDDO                                                      40.80
            ENDIF                                                         40.80
          ENDIF
!
  100     IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) EXIT                        40.31
          CUOPS => CUOPS%NEXTOPS                                          40.31
        END DO                                                            40.31
      ENDIF
!
      RETURN
!   * end of subroutine SPRCON *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWRBC ( COMPDA )                                         40.31 30.90
!                                                                      *
!***********************************************************************
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_GENARR
      USE M_PARALL
      USE SwanGriddata                                                    40.80

      REAL    COMPDA(MCGRD,MCMVAR)                                        30.21

      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SWRBC')
!
      IF (ITEST .GE. 100 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,*) '  ', ICUR, IGTYPE(2)
        WRITE(PRINTF,*) '******** In subroutine SWRBC *******'            30.21
        IF (ICUR .EQ. 1 ) THEN
          WRITE(PRINTF,51)
        ELSE
          WRITE(PRINTF,50)
        ENDIF
      ENDIF
 50   FORMAT('P.index',5X,'coord.',14X,'depth')
 51   FORMAT('P.index',5X,'coord.',14X,'depth',13X,'UX',13X,'UY')

 !
!     *** The arrays start to be filled in the second value ***
!     *** because in COMPDA(1,"variable"), is the default   ***
!     *** value for land points    version 30.21            ***
!
!     ***  Default values for land point  ***
      COMPDA(1,JDP2) = -1.                                                40.00
      IF (JDP1.GT.1) COMPDA(1,JDP1) = -1.                                 40.00
      IF (JDP3.GT.1) COMPDA(1,JDP3) = -1.                                 40.00
      IF (VARWLV) THEN
        COMPDA(1,JWLV2) = 0.                                              40.00
!       next two lines only for nonstat water level
        IF (JWLV1.GT.1) COMPDA(1,JWLV1) = 0.                              40.00
        IF (JWLV3.GT.1) COMPDA(1,JWLV3) = 0.                              40.00
      ENDIF
      IF (ICUR.GT.0) THEN
        COMPDA(1,JVX2) = 0.                                               40.00
        COMPDA(1,JVY2) = 0.                                               40.00
        IF (JVX1.GT.1) COMPDA(1,JVX1) = 0.                                40.00
        IF (JVY1.GT.1) COMPDA(1,JVY1) = 0.                                40.00
        IF (JVX3.GT.1) COMPDA(1,JVX3) = 0.                                40.00
        IF (JVY3.GT.1) COMPDA(1,JVY3) = 0.                                40.00
      ENDIF
      IF (VARFR) THEN
        COMPDA(1,JFRC2) = 0.                                              40.00
        COMPDA(1,JFRC3) = 0.                                              40.00
      ENDIF
      IF (VARWI) THEN
        COMPDA(1,JWX2) = 0.                                               40.00
        COMPDA(1,JWY2) = 0.                                               40.00
        IF (JWX3.GT.1) COMPDA(1,JWX3) = 0.                                40.00
        IF (JWY3.GT.1) COMPDA(1,JWY3) = 0.                                40.00
      ENDIF
      IF (VARAST) THEN
        COMPDA(1,JASTD2) = 0.                                             40.03
        COMPDA(1,JASTD3) = 0.                                             40.03
      ENDIF
      IF (VARNPL) THEN                                                    40.55
        COMPDA(1,JNPLA2) = 0.                                             40.55
        COMPDA(1,JNPLA3) = 0.                                             40.55
      ENDIF
      COMPDA(1,JBOTLV) = 0.                                               40.65
!
!     --- structured grid                                                 40.80
!
      DO 20 IX = 1, MXC
        DO 10 IY = 1, MYC
          INDX = KGRPNT(IX,IY)
          IF (INDX.LE.0 .OR. INDX.GT.MCGRD) THEN                          30.60
            CALL MSGERR (3, 'Grid error in subr. SWRBC')                  30.60
            WRITE (PRINTF, 8) IX, IY, INDX, MCGRD                         30.60
   8        FORMAT (' IX, IY, INDX, MCGRD: ', 4I7)                        30.60
            GOTO 10                                                       30.60
          ENDIF                                                           30.60
          IF (INDX.EQ.1) GOTO 10                                          30.60
!
          XP = XCGRID(IX,IY)                                              30.72
          YP = YCGRID(IX,IY)                                              30.72
!
!         ***** compute depth and water level *****
          DEP = SVALQI (XP, YP, 1, DEPTH, 1, IX, IY)                      40.31 30.90
          COMPDA(INDX,JBOTLV) = DEP                                       40.65
!
          IF (VARWLV) THEN                                                20.38
            WLVL = SVALQI (XP, YP, 7, WLEVL, 1 ,IX ,IY)                   40.31 30.90
            COMPDA(INDX,JWLV2) = WLVL
            IF (JWLV1.GT.1) COMPDA(INDX,JWLV1) = WLVL                     40.00
            IF (JWLV3.GT.1) COMPDA(INDX,JWLV3) = WLVL                     40.00
            DEP = DEP + WLVL
          ENDIF
!         add constant water level
          DEPW = DEP + WLEV                                               30.70
          COMPDA(INDX,JDP2) = DEPW
!         ***In this step the water level at T+DT is copied to the  ***
!         ***water level at T (Only for first time computation)     ***
          IF (JDP1.GT.1) COMPDA(INDX,JDP1) = DEPW                         40.00
          IF (JDP3.GT.1) COMPDA(INDX,JDP3) = DEPW                         40.00
!
!         ***** compute current velocity *****
!
          IF (ICUR.EQ.1 .AND. IGTYPE(2) .GE. 1) THEN                      30.82
            IF (DEPW.GT.0.) THEN
              UU  = SVALQI (XP, YP, 2, UXB, 0 ,IX ,IY)
              VV  = SVALQI (XP, YP, 3, UYB, 0 ,IX ,IY)
              VTOT = SQRT (UU*UU + VV*VV)
              CGMAX = PNUMS(18)*SQRT(GRAV*DEPW)
              IF (VTOT .GT. CGMAX) THEN
                CGFACT = CGMAX / VTOT
                UU = UU * CGFACT
                VV = VV * CGFACT
                IF (ERRPTS.GT.0.AND.IAMMASTER) THEN                       40.95 40.30
                  WRITE (ERRPTS, 211) IX+MXF-1, IY+MYF-1, 1
 211              FORMAT (I4, 1X, I4, 1X, I2)
                ENDIF
              ENDIF
              COMPDA(INDX,JVX2) =  UU*COSVC + VV*SINVC
              COMPDA(INDX,JVY2) = -UU*SINVC + VV*COSVC
            ELSE
              COMPDA(INDX,JVX2) =  0.
              COMPDA(INDX,JVY2) =  0.
            ENDIF
            IF (JVX1.GT.1) COMPDA(INDX,JVX1) = COMPDA(INDX,JVX2)          40.00
            IF (JVY1.GT.1) COMPDA(INDX,JVY1) = COMPDA(INDX,JVY2)          40.00
            IF (JVX3.GT.1) COMPDA(INDX,JVX3) = COMPDA(INDX,JVX2)          40.00
            IF (JVY3.GT.1) COMPDA(INDX,JVY3) = COMPDA(INDX,JVY2)          40.00
          ENDIF

          IF (ITEST .GE. 100 .OR. INTES .GE. 30) THEN
            IF (ICUR .EQ. 1 ) THEN
              WRITE(PRINTF,31)KGRPNT(IX,IY),XP,YP,COMPDA(INDX,JDP2),
     &        COMPDA(INDX,JVX2),COMPDA(INDX,JVY2)
            ELSE
              WRITE(PRINTF,30)KGRPNT(IX,IY),XP,YP,COMPDA(INDX,JDP2)
            ENDIF
          ENDIF
 30       FORMAT(1X,I5,2X,3(E11.4,2X))
 31       FORMAT(1X,I5,2X,5(E11.4,2X))
!
!         ***** compute variable friction coefficient *****
!
          IF (VARFR) THEN
             FRI = SVALQI (XP, YP, 4, FRIC, 1 ,IX ,IY)                    40.31 30.90
             COMPDA(INDX,JFRC2) = FRI                                     40.00
             COMPDA(INDX,JFRC3) = FRI                                     40.00
          ENDIF
!
!         ***** compute variable wind velocity *****
!
          IF (VARWI) THEN
             UU  = SVALQI (XP, YP, 5, WXI, 0 ,IX ,IY)                     40.31 30.90
             VV  = SVALQI (XP, YP, 6, WYI, 0 ,IX ,IY)                     40.31 30.90
             COMPDA(INDX,JWX2) =  UU*COSWC + VV*SINWC
             COMPDA(INDX,JWY2) = -UU*SINWC + VV*COSWC
            IF (JWX3.GT.1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          40.00
            IF (JWY3.GT.1) COMPDA(INDX,JWY3) = COMPDA(INDX,JWY2)          40.00
          ENDIF
!
!     ***** compute variable air-sea temperature difference *****         40.03
!
          IF (VARAST) THEN
             ASTD = SVALQI (XP, YP, 10, ASTDF, 1 ,IX ,IY)                 40.31 40.03
             COMPDA(INDX,JASTD2) = ASTD                                   40.03
             COMPDA(INDX,JASTD3) = ASTD                                   40.03
          ENDIF
!
!     ***** compute number of plants per square meter *****               40.55
!
          IF (VARNPL) THEN
             XNPL = SVALQI (XP, YP, 11, NPLAF, 1 ,IX ,IY)                 40.55
             COMPDA(INDX,JNPLA2) = XNPL                                   40.55
             COMPDA(INDX,JNPLA3) = XNPL                                   40.55
          ENDIF
!
   10   CONTINUE
   20 CONTINUE
!
!     --- unstructured grid                                               40.80
!
      DO INDX = 1, nverts
!
         XP = xcugrd(INDX)                                                40.80
         YP = ycugrd(INDX)                                                40.80
!
!        ***** compute depth and water level *****
!
         IF ( IGTYPE(1).EQ.3 ) THEN
            DEP = DEPTH(INDX)
         ELSE
            DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)                      40.80
         ENDIF
         COMPDA(INDX,JBOTLV) = DEP                                        40.80
!
         IF (VARWLV) THEN                                                 40.80
            IF ( IGTYPE(7).EQ.3 ) THEN
               WLVL = WLEVL(INDX)
            ELSE
               WLVL = SVALQI (XP, YP, 7, WLEVL, 1, 0, 0)                  40.80
            ENDIF
            COMPDA(INDX,JWLV2) = WLVL
            IF (JWLV1.GT.1) COMPDA(INDX,JWLV1) = WLVL                     40.80
            IF (JWLV3.GT.1) COMPDA(INDX,JWLV3) = WLVL                     40.80
            DEP = DEP + WLVL
         ENDIF
!        add constant water level
         DEPW = DEP + WLEV                                                40.80
         COMPDA(INDX,JDP2) = DEPW
!        ***In this step the water level at T+DT is copied to the  ***
!        ***water level at T (Only for first time computation)     ***
         IF (JDP1.GT.1) COMPDA(INDX,JDP1) = DEPW                          40.80
         IF (JDP3.GT.1) COMPDA(INDX,JDP3) = DEPW                          40.80
!
!        ***** compute current velocity *****
!
         IF (ICUR.EQ.1 .AND. IGTYPE(2) .GE. 1) THEN                       40.80
            IF (DEPW.GT.0.) THEN
              IF ( IGTYPE(2).EQ.3 ) THEN
                 UU = UXB(INDX)
              ELSE
                 UU = SVALQI (XP, YP, 2, UXB, 0, 0, 0)                    40.80
              ENDIF
              IF ( IGTYPE(3).EQ.3 ) THEN
                 VV = UYB(INDX)
              ELSE
                 VV = SVALQI (XP, YP, 3, UYB, 0, 0, 0)                    40.80
              ENDIF
              VTOT = SQRT (UU*UU + VV*VV)
              CGMAX = PNUMS(18)*SQRT(GRAV*DEPW)
              IF (VTOT .GT. CGMAX) THEN
                CGFACT = CGMAX / VTOT
                UU = UU * CGFACT
                VV = VV * CGFACT
                IF (ERRPTS.GT.0) THEN
                   WRITE (ERRPTS, 212) INDX, 1
 212               FORMAT (I4, 1X, I2)
                ENDIF
              ENDIF
              COMPDA(INDX,JVX2) =  UU*COSVC + VV*SINVC
              COMPDA(INDX,JVY2) = -UU*SINVC + VV*COSVC
            ELSE
              COMPDA(INDX,JVX2) =  0.
              COMPDA(INDX,JVY2) =  0.
            ENDIF
            IF (JVX1.GT.1) COMPDA(INDX,JVX1) = COMPDA(INDX,JVX2)          40.80
            IF (JVY1.GT.1) COMPDA(INDX,JVY1) = COMPDA(INDX,JVY2)          40.80
            IF (JVX3.GT.1) COMPDA(INDX,JVX3) = COMPDA(INDX,JVX2)          40.80
            IF (JVY3.GT.1) COMPDA(INDX,JVY3) = COMPDA(INDX,JVY2)          40.80
         ENDIF
!
!         ***** compute variable friction coefficient *****
!
         IF (VARFR) THEN
            IF ( IGTYPE(4).EQ.3 ) THEN
               FRI = FRIC(INDX)
            ELSE
               FRI = SVALQI (XP, YP, 4, FRIC, 1, 0, 0)                    40.80
            ENDIF
            COMPDA(INDX,JFRC2) = FRI                                      40.80
            COMPDA(INDX,JFRC3) = FRI                                      40.80
         ENDIF
!
!        ***** compute variable wind velocity *****
!
         IF (VARWI) THEN
            IF ( IGTYPE(5).EQ.3 ) THEN
               UU = WXI(INDX)
            ELSE
               UU = SVALQI (XP, YP, 5, WXI, 0, 0, 0)                      40.80
            ENDIF
            IF ( IGTYPE(6).EQ.3 ) THEN
               VV = WYI(INDX)
            ELSE
               VV = SVALQI (XP, YP, 6, WYI, 0, 0, 0)                      40.80
            ENDIF
            COMPDA(INDX,JWX2) =  UU*COSWC + VV*SINWC
            COMPDA(INDX,JWY2) = -UU*SINWC + VV*COSWC
            IF (JWX3.GT.1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          40.80
            IF (JWY3.GT.1) COMPDA(INDX,JWY3) = COMPDA(INDX,JWY2)          40.80
         ENDIF
!
!     ***** compute variable air-sea temperature difference *****         40.80
!
         IF (VARAST) THEN
            IF ( IGTYPE(10).EQ.3 ) THEN
               ASTD = ASTDF(INDX)
            ELSE
               ASTD = SVALQI (XP, YP, 10, ASTDF, 1, 0, 0)                 40.80
            ENDIF
            COMPDA(INDX,JASTD2) = ASTD                                    40.80
            COMPDA(INDX,JASTD3) = ASTD                                    40.80
         ENDIF
!
!     ***** compute number of plants per square meter *****               40.80
!
          IF (VARNPL) THEN
             IF ( IGTYPE(11).EQ.3 ) THEN
                XNPL = NPLAF(INDX)
             ELSE
                XNPL = SVALQI (XP, YP, 11, NPLAF, 1 ,0, 0)                40.80
             ENDIF
             COMPDA(INDX,JNPLA2) = XNPL                                   40.80
             COMPDA(INDX,JNPLA3) = XNPL                                   40.80
          ENDIF
!
      ENDDO                                                               40.80
!
!     *** initialise setup and saved depth ***                            32.02
!
      IF (LSETUP.GT.0) THEN                                               32.02
        DO INDX = 1, MCGRD                                                32.02
          COMPDA(INDX,JSETUP) =  0.                                       32.02
          COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)                         32.02
        ENDDO                                                             32.02
      ENDIF                                                               32.02
!
!     --- initialize HSIBC                                                40.41
      COMPDA(:,JHSIBC) = 0.                                               40.41
!
!     --- initialize ZELEN and USTAR                                      40.41
      COMPDA(:,JZEL  ) = 2.E-33                                           40.41
      COMPDA(:,JUSTAR) = 1.E-15                                           40.41
!
!     --- initialize UBOT and TMBOT                                       40.94
      COMPDA(:,JUBOT) = 0.                                                40.94
      IF (JPBOT.GT.1) COMPDA(:,JPBOT) = 0.                                40.94
!
      RETURN
! * end of subroutine SWRBC *
      END
!***********************************************************************
!                                                                      *
      REAL FUNCTION SVALQI (XP, YP, IGRID, ARRINP, ZERO ,IXC ,IYC)        30.21
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE M_PARALL

      IMPLICIT NONE
!
      INTEGER  IGRID, IXC, IYC, ZERO
      REAL     ARRINP(*), XP, YP
      INTEGER  IB1, IENT, II, JB1, IXCGL, IYCGL
      REAL     IXB, IYB, SXB1, SXB2, SYB1, SYB2
      REAL     SUMWEXC, SUMWREG, WF1, WF2, WF3, WF4                       40.04
      LOGICAL  INGRD
      LOGICAL  EQREAL
      SAVE IENT                                                           30.72
      DATA IENT/0/
      CALL STRACE (IENT, 'SVALQI')

!     --- take global indices instead of local ones                       40.30

      IXCGL = IXC + MXF - 1                                               40.30
      IYCGL = IYC + MYF - 1                                               40.30
!
!     ***    Two different procedures in funcion of       ***
!     ***    grid type: regular or curvilinear (staggered)*** ver 30.21
!
      IF (IGTYPE(IGRID) .EQ. 1) THEN                                      30.21
!
!     Regular grid:
!
        IXB = ( (XP-XPG(IGRID))*COSPG(IGRID) +
     &          (YP-YPG(IGRID))*SINPG(IGRID) ) / DXG(IGRID)
!
        INGRD = .TRUE.
        IF (IXB .LE. 0.) THEN                                             30.60
          IB1   = 1
          SXB2  = 0.
          IF (IXB.LT.-0.1) INGRD = .FALSE.                                20.3x
        ELSE IF (IXB .GE. FLOAT(MXG(IGRID)-1)) THEN                       30.60
          IB1   = MXG(IGRID)-1
          SXB2  = 1.
          IF (IXB.GT.FLOAT(MXG(IGRID))-0.9) INGRD = .FALSE.               20.3x
        ELSE
          IB1   = INT(IXB)
          SXB2  = IXB-REAL(IB1)
          IB1   = IB1+1
        ENDIF
        IF (MYG(IGRID).GT.1) THEN                                         32.03
          IYB = (-(XP-XPG(IGRID))*SINPG(IGRID) +                          30.82
     &            (YP-YPG(IGRID))*COSPG(IGRID) ) / DYG(IGRID)             30.82
          IF (IYB .LE. 0.) THEN                                           30.60
            JB1   = 1
            SYB2  = 0.
            IF (IYB.LT.-0.1) INGRD = .FALSE.
          ELSE IF (IYB .GE. FLOAT(MYG(IGRID)-1)) THEN                     30.60
            JB1   = MYG(IGRID)-1
            SYB2  = 1.
            IF (IYB.GT.FLOAT(MYG(IGRID))-0.9) INGRD = .FALSE.
          ELSE
            JB1   = INT(IYB)
            SYB2  = IYB-REAL(JB1)
            JB1   = JB1+1
          ENDIF
        ENDIF                                                             32.03
!
!       evaluate SVALQI (2D-mode):
!
        IF (.NOT.INGRD .AND. ZERO.EQ.0) THEN
          SVALQI = 0.
        ELSE IF (MYG(IGRID).GT.1) THEN                                    32.03
          SXB1   = 1.- SXB2
          SYB1   = 1.- SYB2
          II     = IB1 + (JB1-1) * MXG(IGRID)
          WF1 = SXB1*SYB1                                                 40.04
          WF2 = SXB1*SYB2                                                 40.04
          WF3 = SXB2*SYB1                                                 40.04
          WF4 = SXB2*SYB2                                                 40.04
          SUMWEXC = 0.                                                    40.04
          IF  (EQREAL(ARRINP(II             ),EXCFLD(IGRID))) THEN        40.04
            SUMWEXC = SUMWEXC + WF1                                       40.04
            WF1 =0.                                                       40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+  MXG(IGRID)),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF2                                       40.04
            WF2=0.                                                        40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+1           ),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF3                                       40.04
            WF3=0.                                                        40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+1+MXG(IGRID)),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF4                                       40.04
            WF4=0.                                                        40.04
          ENDIF                                                           40.04
          SUMWREG = 1. -SUMWEXC                                           40.04
!
          IF (SUMWEXC.GE.SUMWREG)   THEN                                  40.04
            SVALQI = EXCFLD(IGRID)                                        40.04
          ELSE                                                            40.04
            SVALQI = ( WF1*ARRINP(II)   + WF2*ARRINP(II+MXG(IGRID))       40.04
     &               + WF3*ARRINP(II+1) + WF4*ARRINP(II+1+MXG(IGRID)) )   40.04
     &               / SUMWREG                                            40.04
          END IF                                                          40.04
        ELSE
!
!       evaluate SVALQI (1D-mode):
!
          SXB1 = 1. - SXB2                                                32.03
          IF (EQREAL(ARRINP(IB1  ),EXCFLD(IGRID)).OR.                     30.82
     &        EQREAL(ARRINP(IB1+1),EXCFLD(IGRID))    ) THEN               30.82
!
!           One of the cornerpoints contains an exception value thus:     30.82
!
            SVALQI = EXCFLD(IGRID)                                        30.82
          ELSE                                                            30.82
            SVALQI = SXB1*ARRINP(IB1)                                     32.03
     &             + SXB2*ARRINP(IB1+1)                                   32.03
          ENDIF                                                           30.82
        ENDIF
      ELSEIF ( IGTYPE(IGRID).EQ.3 ) THEN                                  40.80
!
!     unstructured grid
!
        CALL SwanInterpolatePoint(SVALQI, XP, YP, ARRINP, EXCFLD(IGRID))
!
      ELSE IF (ABS(STAGX(IGRID)) .LT. 0.01 .AND.                          32.03
     &         ABS(STAGY(IGRID)) .LT. 0.01) THEN                          32.03
!
!     Curvi-linear and non-staggered input grid:                          32.03
!
        IB1   = IXCGL
        JB1   = IYCGL
        II     = IB1 + (JB1-1) * MXG(IGRID)
        SVALQI = ARRINP(II)
      ELSE
!
!     Curvi-linear and staggered input grid:                              32.03
!
        INGRD = .TRUE.
        IF (IXCGL .EQ. 1) THEN
          IB1   = 1
          SXB2  = 0.
          IF (STAGY(IGRID) .GT. 0.) INGRD = .FALSE.
        ELSE IF (IXCGL .GT. MXG(IGRID)-1) THEN
          IB1   = MXG(IGRID)-1
          SXB2  = 1.
          IF (STAGY(IGRID) .GT. 0.) INGRD = .FALSE.
        ELSE
          IB1   = IXCGL + 1
          SXB2  = 1. - STAGX(IGRID)
        ENDIF
        IF (IYCGL .EQ. 1) THEN
          JB1   = 1
          SYB2  = 0.
          IF (STAGX(IGRID) .GT. 0.) INGRD = .FALSE.
        ELSE IF (IYCGL .GT. MYG(IGRID)-1) THEN
          JB1   = MYG(IGRID)-1
          SYB2  = 1.
          IF (STAGY(IGRID) .GT. 0.) INGRD = .FALSE.
        ELSE
          JB1   = IYCGL + 1
          SYB2  =1. - STAGY(IGRID)
        ENDIF
!
!       evaluate SVALQI (2D-mode):
!
        IF (.NOT.INGRD .AND. ZERO.EQ.0) THEN
          SVALQI = 0.
        ELSE
          SXB1   = STAGX(IGRID)
          SYB1   = STAGY(IGRID)
          II     = IB1 + (JB1-1) * MXG(IGRID)
          WF1 = SXB1*SYB1                                                 40.04
          WF2 = SXB1*SYB2                                                 40.04
          WF3 = SXB2*SYB1                                                 40.04
          WF4 = SXB2*SYB2                                                 40.04
          SUMWEXC = 0.                                                    40.04
          IF (EQREAL(ARRINP(II             ),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF1                                       40.04
            WF1 =0.                                                       40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+  MXG(IGRID)),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF2                                       40.04
            WF2=0.                                                        40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+1           ),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF3                                       40.04
            WF3=0.                                                        40.04
          ENDIF                                                           40.04
          IF (EQREAL(ARRINP(II+1+MXG(IGRID)),EXCFLD(IGRID))) THEN         40.04
            SUMWEXC = SUMWEXC + WF4                                       40.04
            WF4=0.                                                        40.04
          ENDIF                                                           40.04
          SUMWREG = 1. -SUMWEXC                                           40.04
!
          IF (SUMWEXC.GE.SUMWREG)   THEN                                  40.04
            SVALQI = EXCFLD(IGRID)                                        40.04
          ELSE                                                            40.04
            SVALQI = ( WF1*ARRINP(II)   + WF2*ARRINP(II+MXG(IGRID))       40.04
     &               + WF3*ARRINP(II+1) + WF4*ARRINP(II+1+MXG(IGRID)) )   40.04
     &               / SUMWREG                                            40.04
          END IF                                                          40.04
        ENDIF
      ENDIF
!
!     ***** test *****
      IF (ITEST .GE. 280)
!     &   WRITE(PRINTF, 6010) SVALQI,IGRID,XP,YP,IXB,IYB,II,ARRINP(II)
! 6010 FORMAT(' SVALQI  IGRID       XP      YP        IXB',
!     &       '       IYB  II  ARRINP(II)', /
!     &      ,E10.3,I3,1X,4E10.3,I4,E10.3)
     &   WRITE(PRINTF, 6010) XP, YP, IXB, IYB, SVALQI
 6010 FORMAT(' Test SVALQI:',5F10.3)
!
      RETURN
!     end of subroutine SVALQI
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SINUPT (PSNAME, XP, YP, XCGRID, YCGRID, KGRPNT, KGRBND)  40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3
!
      INTEGER KGRPNT(MXC,MYC), KGRBND(*)
      REAL    XP, YP
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)
      CHARACTER PSNAME *(*)
      LOGICAL SINBTG, SINCMP
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINUPT')
!
      IF (.NOT. SINBTG (XP,YP) ) THEN
         CALL MSGERR(1,'(corner)point outside bottom grid')
         WRITE (PRINTF, 6010) PSNAME, XP+XOFFS, YP+YOFFS
      ENDIF
      IF (.NOT.SINCMP (XP, YP, XCGRID, YCGRID, KGRPNT, KGRBND)) THEN      40.00
        CALL MSGERR(1,'(corner)point outside comp. grid')
        WRITE (PRINTF, 6010) PSNAME, XP+XOFFS, YP+YOFFS
      ENDIF
 6010 FORMAT('       Set of output locations: ',A8,
     &       '  coordinates:', 2F10.2)
!
      RETURN
!     end of subroutine SINUPT *
      END
!
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION SINBTG (XP, YP)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM2                                                         40.41

      REAL    XB, YB, XLENB, YLENB, BTOL
      DATA  IENT /0/
      CALL  STRACE (IENT,'SINBTG')
!
      SINBTG = .TRUE.
      IF ( IGTYPE(1).NE.1 ) RETURN
!
      XLENB = (MXG(1)-1)*DXG(1)
      YLENB = (MYG(1)-1)*DYG(1)
      BTOL  = 0.01 * (XLENB+YLENB)
!
!     ***** compute bottom grid coordinates from problem coordinates ****
!
      XB =  (XP-XPG(1))*COSPG(1) + (YP-YPG(1))*SINPG(1)
      YB = -(XP-XPG(1))*SINPG(1) + (YP-YPG(1))*COSPG(1)
!
!     ***** check location of point *****
      IF (XB .LT. -BTOL) SINBTG = .FALSE.
      IF (XB .GT. XLENB+BTOL) SINBTG = .FALSE.
      IF (YB .LT. -BTOL) SINBTG = .FALSE.
      IF (YB .GT. YLENB+BTOL) SINBTG = .FALSE.
!
      RETURN
!   * end of subroutine SINBTG *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION SINCMP (XP, YP ,XCGRID ,YCGRID ,KGRPNT, KGRBND)    40.00
!                                                                      *
!***********************************************************************
!
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
!PUN      USE SIZES, ONLY: MNPROC
!
      INTEGER KGRPNT(MXC,MYC), KGRBND(*)                                  40.00

      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
      REAL    XP,     YP
      REAL    CTOL

      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINCMP')
!
!     *** Different procedure depending on grid type **
      IF (OPTG .EQ. 1) THEN                                               30.21
!
!       regular grid: compute comp. coordinates from problem coordinates
!
        XC   =  (XP-XPC)*COSPC+(YP-YPC)*SINPC
        YC   = -(XP-XPC)*SINPC+(YP-YPC)*COSPC
!       XC and YC are in m
!
!       ***** check for location *****
        SINCMP = .TRUE.
        CTOL   = 0.01 * (XCLEN+YCLEN)                                     40.00
        IF (XC .LT. -CTOL) SINCMP = .FALSE.                               40.00
        IF (XC .GT. XCLEN+CTOL) SINCMP = .FALSE.                          40.00
        IF (YC .LT. -CTOL) SINCMP = .FALSE.                               40.00
        IF (YC .GT. YCLEN+CTOL) SINCMP = .FALSE.                          40.00
      ELSE IF (OPTG .EQ. 3) THEN                                          40.80
!
!       curvilinear grid
!
        CALL CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID, KGRBND)      40.00
!       XC and YC are nondimensional; equivalent to grid index
!
!       ***** check for location *****
        SINCMP = .TRUE.
        IF (XC .LT. -0.01) SINCMP = .FALSE.                               40.00
        IF (XC .GT. REAL(MXC-1)+0.01) SINCMP = .FALSE.                    40.00
        IF (YC .LT. -0.01) SINCMP = .FALSE.                               40.00
        IF (YC .GT. REAL(MYC-1)+0.01) SINCMP = .FALSE.                    40.00
      ELSE IF (OPTG.EQ.5) THEN                                            40.80
!
!       unstructured grid
!
        SINCMP = .TRUE.                                                   40.80
        CALL SwanFindPoint ( XP, YP, K )                                  40.80
        IF ( K.LT.0 ) SINCMP = .FALSE.                                    40.80
!
!PUN!       --- check if output location is in global subdomain               40.95
!PUN!
!PUN        IF ( MNPROC.GT.1 .AND. .NOT.SINCMP ) THEN                         40.95
!PUN           IF ( XP.GE.XCGMIN .AND. XP.LE.XCGMAX .AND.                     40.95
!PUN     &          YP.GE.YCGMIN .AND. YP.LE.YCGMAX ) SINCMP = .TRUE.         40.95
!PUN        ENDIF                                                             40.95
!PUN!
      ENDIF
!
!     --- check if output location is in global subdomain                 40.31
!
      IF ( PARLL .AND. .NOT.SINCMP ) THEN                                 40.31
         IF ( XP.GE.XCGMIN .AND. XP.LE.XCGMAX .AND.                       40.31
     &        YP.GE.YCGMIN .AND. YP.LE.YCGMAX ) SINCMP = .TRUE.           40.31
      END IF                                                              40.31
!
      RETURN
!   * end of subroutine SINCMP *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE WRTEST (NAME, NA, IARR, RARR)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      INTEGER   IARR(*), NA                                               30.72
      REAL      RARR(*)
      CHARACTER NAME *(*)
!
      WRITE (PRINTF, 10) NAME, (IARR(II), II=1,NA)
  10  FORMAT (1X, A, 10(1X, I8))
      WRITE (PRINTF, 20) (RARR(II), II=1,NA)
  20  FORMAT (10(1X, E12.4))
      RETURN
! * end of subroutine WRTEST *
      END
!********************************************************************
!
      SUBROUTINE ERRCHK                                                   40.31 40.00
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_GENARR                                                        40.41
      CHARACTER*80 MSGSTR
      LOGICAL  EQREAL                                                     40.53
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'ERRCHK')

!     --- choose scheme for stationary or nonstationary computation

      IF (NSTATC.EQ.1) THEN                                               40.31
        PROPSC = PROPSN                                                   40.31
      ELSE                                                                40.31
        PROPSC = PROPSS                                                   40.31
      ENDIF                                                               40.31
!
!     --- set default setting for stopping criterion                      40.80
!
      IF ( PNUMS(21).EQ.-9. ) THEN                                        40.80
         IF (OPTG.NE.5) THEN                                              40.80
            PNUMS(21) = 0.                                                40.80
         ELSE                                                             40.80
            PNUMS(21) = 1.                                                40.80
            PNUMS( 1) = 0.01                                              40.80
            PNUMS( 2) = 0.005                                             40.80
            PNUMS( 4) = 99.5                                              40.80
            PNUMS(15) = 0.005                                             40.80
         ENDIF                                                            40.80
      ENDIF                                                               40.80
!
!     --- reset PWIND(9) and PWIND(17) if RHO has a user value            40.96
!
      PWIND( 9) = PWIND(16)/RHO                                           40.96
      PWIND(17) = RHO                                                     40.96
!
!     -----------------------------------------------------------------
!
!     *** WARNINGS AND ERROR MESSAGES ***
!
!     -----------------------------------------------------------------
!
!     *** Check formulation for whitecapping ***
!
      IF ( IWCAP.GT.7 ) THEN
         WRITE (MSGSTR, '(A,I2,A)')
     &               'Unknown method for whitecapping (IWCAP=',IWCAP,')'
         CALL MSGERR( 1, TRIM(MSGSTR) )
         CALL MSGERR( 1,
     &     'Whitecapping according to Komen et al. (1984) will be used')
         IWCAP     = 1
         PWCAP(1)  = 2.36E-5
         PWCAP(2)  = 3.02E-3
         PWCAP(9)  = 2.
         PWCAP(10) = 0.
         PWCAP(11) = 1.
      ENDIF
!
!     *** WAM cycle 3 physics ***
!
      IF  ( (IWIND .EQ. 3 .OR. IWIND .EQ. 5) .AND. IWCAP .NE. 1
     &      .AND. IWCAP.NE.7                                              40.53
     &    ) THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Komen et al. (1984) for wind option G3/YAN')       20.74
      ENDIF

      IF  ( IWCAP.EQ.7 .AND. IWIND.NE.5 ) THEN                            40.53
        CALL MSGERR(1,'Activate wind option Yan (1987) in case of')       40.53
        CALL MSGERR(1,'Alves & Banner (2003) white-capping method')       40.53
      END IF                                                              40.53
!
!     *** WAM cycle 4 physics ***
!
      IF  ( IWIND .EQ. 4 .AND. IWCAP .NE. 2) THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Janssen (1991) for wind option JANS        ')
      ENDIF
!
!     *** check numerical scheme in presence of a current ***
!
      IF ( ICUR .EQ. 1 ) THEN
        IF ( PNUMS(6) .EQ. 0. ) THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in theta space ')
          CALL MSGERR(1,'-> set CDD = 1.')
          WRITE(PRINTF,*)
        ENDIF
        IF ( PNUMS(7) .EQ. 0. ) THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in sigma space ')
          CALL MSGERR(1,'-> set CSS = 1.')
          WRITE(PRINTF,*)
        ENDIF
      END IF
!
!     check combination of REPeating option and grid type and dimension   33.09
!
      IF (KREPTX.GT.0) THEN                                               33.08
!       --- "GE" changed to "EQ" since OPTG.EQ.3 FOR CURVILINEAR          40.08
        IF (OPTG.EQ.3)                                                    40.08
     &  CALL MSGERR (3, 'Curvilinear grid cannot be REPeating')
        IF (OPTG.EQ.5)                                                    40.80
     &  CALL MSGERR (3, 'Unstructured grid cannot be REPeating')          40.80
        IF (PROPSC.EQ.1 .AND. MXC.LT.1)
     &  CALL MSGERR (3, 'MXC must be >=1 for REPeating option')
        IF (PROPSC.EQ.2 .AND. MXC.LT.2)                                   33.10
     &  CALL MSGERR (3, 'MXC must be >=2 for REPeating option')           33.10
        IF (PROPSC.EQ.3 .AND. MXC.LT.3)                                   33.09
     &  CALL MSGERR (3, 'MXC must be >=3 for REPeating option')           33.09
      ENDIF
!
      IF (PROPSC.EQ.2 .AND. NSTATC.GT.0) THEN                             33.10
        CALL MSGERR (3, 'SORDUP scheme only in stationary run')           33.10
      ENDIF
      IF (PROPSC.EQ.3 .AND. NSTATC.EQ.0) THEN                             33.08
        CALL MSGERR (3, 'S&L scheme not in stationary run')               33.08
      ENDIF
!
!     --- A warning about curvilinear and S&L scheme                      40.08
      IF ((PROPSC.EQ.3).AND.(OPTG.EQ.3)) THEN                             40.08
         CALL MSGERR(1,'the S&L scheme (higher order nonstationary')      40.08
         CALL MSGERR(1,'IS NOT fully implemented for curvilinear')        40.08
         CALL MSGERR(1,'coordinates. This may or may not be noticeable')  40.08
         CALL MSGERR(1,'in simulations. DX and DY are approximated')      40.08
         CALL MSGERR(1,'with DX and DY of two nearest cells.')            40.08
         CALL MSGERR(1,'Note that SORDUP (higher order stationary)')      40.08
         CALL MSGERR(1,'IS fully implemented for curvilinear coord.')     40.08
         CALL MSGERR(1,'and differences are usually negligible.')         40.08
      ENDIF                                                               40.08
!
!     Here the various problems with quadruplets are checked
!     The combination of quadruplets and sectors is an error
!     in the calculation of quadruplets when the SECTOR option is
!     used in the CGRID command. This error should be corrected
!     in the future
!
      IF (IWIND.EQ.3 .OR. IWIND.EQ.4) THEN                                30.60
        IF (IQUAD .EQ. 0) THEN
          CALL MSGERR(2,'Quadruplets should be activated when SWAN  ')    30.60
          CALL MSGERR(2,'is running in a third generation mode and  ')    30.60
          CALL MSGERR(2,'wind is present                            ')    30.60
        ENDIF
      ENDIF
!
      IF (IQUAD .GE. 1) THEN                                              30.72
!
       IF (.NOT. FULCIR) THEN                                             30.72
        IF ((SPDIR2-SPDIR1) .LT. (PI/12.)) THEN                           30.72
          CALL MSGERR(2,'A combination of using quadruplets with a'    )  30.72
          CALL MSGERR(2,'sector of less than 30 degrees should be'     )  30.72
          CALL MSGERR(2,'avoided at all times, it is likely to produce')  30.72
          CALL MSGERR(2,'unreliable results and unexpected errors.'    )  30.72
          CALL MSGERR(2,'Refer to the manual (CGRID) for details'      )  30.72
        ELSE                                                              30.72
          CALL MSGERR(1,'It is not recommended to use quadruplets'     )  30.72
          CALL MSGERR(1,'in combination with calculations on a sector.')  30.72
          CALL MSGERR(1,'Refer to the manual (CGRID) for details'      )  30.72
        END IF                                                            30.72
       END IF                                                             30.72
!
       IF (IWIND.EQ.0) THEN                                               30.72
         CALL MSGERR(2,'It is not recommended to use quadruplets'     )   30.72
         CALL MSGERR(2,'in combination with zero wind conditions.'    )   30.72
       END IF                                                             30.72
!
       IF (MSC .EQ. 3 ) THEN
         CALL MSGERR(4,'Do not activate quadruplets for boundary ')
         CALL MSGERR(4,'option BIN -> use other option           ')
         RETURN
       END IF
!
      END IF                                                              30.72
!
!     Check whether limiter should be de-activated
!
      IF (IQUAD.EQ.0 .AND. PNUMS(20).LT.100.) THEN                        40.41
         CALL MSGERR(1,
     &              'Limiter is de-activated in case of no quadruplets')  40.41
         PNUMS(20) = 1.E+20
      END IF
!
!     Check resolution in frequency-space when DIA is used                40.41
!
      IF (IQUAD.GT.0 .AND. IQUAD.LE.3 .OR. IQUAD.EQ.8) THEN               40.41
         GAMMA = EXP(ALOG(SHIG/SLOW)/REAL(MSC-1))                         40.31
         IF (ABS(GAMMA-1.1).GT.0.055) THEN                                40.31
            CALL MSGERR(1,                                                40.31
     &           'relative frequency resolution (df/f) deviates more')    40.31
            CALL MSGERR(1,                                                40.31
     &           'than 5% from 10%-resolution. This may be problematic')  40.31
            CALL MSGERR(1,                                                40.31
     &           'when quadruplets are approximated by means of DIA.')    40.31
         END IF                                                           40.31
      END IF                                                              40.41
!
!     When using triads MSC must be less than 200!                        30.72
!
      IF ((ITRIAD.GT.0).AND.(MSC.GT.200)) THEN                            30.72
         CALL MSGERR(4,'When triads are active the number of     ')       30.72
         CALL MSGERR(4,'directions chosen in the CGRID command   ')       30.72
         CALL MSGERR(4,'must be less than 200                    ')       30.72
         RETURN
      END IF                                                              30.72
!
!     When Alves and Banner and XNL are applied change the parameters     40.53
!
      IF ( IWCAP.EQ.7 .AND. (IQUAD.EQ.51 .OR. IQUAD.EQ.52 .OR.            40.53
     &                       IQUAD.EQ.53) ) THEN                          40.53
         IF (EQREAL(PWCAP( 1), 5.0E-5)) PWCAP( 1) = 5.0E-5                40.53
         IF (EQREAL(PWCAP(12),1.75E-3)) PWCAP(12) = 1.95E-3               40.53
         IF (EQREAL(PWCAP(10),    4.0)) PWCAP(10) = 4.                    40.53
         IF (EQREAL(PWCAP( 9),    0.0)) PWCAP( 9) = 0.                    40.53
         IF (EQREAL(PWCAP(11),    0.0)) PWCAP(11) = 0.                    40.53
      END IF                                                              40.53
!
      IF ( ITEST .GE. 120 ) THEN
        WRITE(PRINTF,3000) IWIND ,IQUAD, ICUR, IWCAP, MSC
3000    FORMAT(' ERRCHK : IWIND QUAD CUR WCAP MSC   : ',5I4)
        IF (IWIND .GT. 0) THEN                                            24/MAR
          DO II = 1, MWIND
            WRITE(PRINTF,30) II,PWIND(II)
 30         FORMAT(' PWIND(',I2,') = ',E11.4)
          ENDDO
        ENDIF
      END IF
!
      RETURN
!     end of subroutine ERRCHK
      END
!*********************************************************************
!                                                                    *
      SUBROUTINE SNEXTI (BSPECS, BGRIDP, COMPDA, AC1   , AC2   ,          40.31
     &                   SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT,          40.31
     &                   XYTST , DEPTH , WLEVL , FRIC  , UXB   ,          40.31
     &                   UYB   , NPLAF , WXI   , WYI   )                  40.55 40.31
!                                                                    *
!*********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC                                                       40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata
!
      INTEGER  BGRIDP(*), XYTST(*), KGRPNT(MXC,MYC)
      REAL     AC1(MDC,MSC,MCGRD)                                         40.00
      REAL     AC2(MDC,MSC,MCGRD)                                         30.21
      REAL     BSPECS(MDC,MSC,NBSPEC,2)                                   40.00
      REAL     COMPDA(MCGRD,MCMVAR)                                       30.72
      REAL     SPCDIR(MDC,6)                                              40.00
      REAL     SPCSIG(MSC)                                                30.01
      REAL     XCGRID(MXC,MYC), YCGRID(MXC,MYC)                           30.21
      REAL     DEPTH(*), WLEVL(*), FRIC(*), UXB(*), UYB(*),               40.31
     &         NPLAF(*),                                                  40.55
     &         WXI(*), WYI(*)                                             40.3
      LOGICAL STPNOW
      INTEGER     IERR
      REAL        TSTVAL(10)                                              40.00
      TYPE(BSPCDAT), POINTER :: CURBFL                                    40.31
      LOGICAL LPB                                                         40.80
      SAVE  IENT
      DATA  IENT/0/
      IF (LTRACE) CALL STRACE(IENT,'SNEXTI')
!
!
!     **   All action densities are shifted from array T+DT
!     **   to the array at time T
!
      IF (NSTATC.EQ.1) THEN                                               30.70
        IF (ITERMX.GT.1
     &       .OR. PROPSC.EQ.3                                             33.09
     &                                 ) THEN                             30.70
          DO 50 IXY = 1, MCGRD                                            30.21
            DO 55 ISS = 1, MSC
              DO 60 IDD = 1, MDC
                AC1(IDD,ISS,IXY) = AC2(IDD,ISS,IXY)                       30.21
 60           CONTINUE
 55         CONTINUE
 50       CONTINUE
        ENDIF                                                             40.00
      ENDIF

!     --- update boundary conditions

      IF (INODE.EQ.MASTER) THEN                                           40.30
         IF (ITEST.GE.80) WRITE (PRTEST,*) ' number of boundary files ',
     &             NBFILS
         CURBFL => FBNDFIL
         DO IBFILE = 1, NBFILS

!          --- read values from boundary file, and interpolate in time

           CALL RBFILE ( SPCSIG, SPCDIR,                                  40.31
     &                   CURBFL%BFILED, CURBFL%BSPLOC,                    40.31
     &                   CURBFL%BSPDIR, CURBFL%BSPFRQ,                    40.31
     &                   BSPECS, XYTST )                                  40.31
           IF (STPNOW()) RETURN                                           34.01
           IF (.NOT.ASSOCIATED(CURBFL%NEXTBSPC)) EXIT                     40.31
           CURBFL => CURBFL%NEXTBSPC                                      40.31

         END DO
      END IF

!     --- scatter array BSPECS to all nodes
      CALL SWBROADC ( BSPECS, 2*MDC*MSC*NBSPEC, SWREAL )                  40.30
      IF (STPNOW()) RETURN                                                40.30

!     --- determine spectra on boundary points of grid

      IF ( NBGRPT.GT.0 ) THEN
         IF (ITEST.GE.80) WRITE(PRTEST,*) ' number of boundary points ',
     &             NBGRPT
         DO IBGRID = 1, NBGRPT
           INDXGR = BGRIDP(6*IBGRID-5)
           IF (BGRIDP(6*IBGRID-4).EQ.1) THEN
!            obtain spectrum in boundary point from interpolation in space
             W1 = 0.001 * REAL(BGRIDP(6*IBGRID-3))
             K1 = BGRIDP(6*IBGRID-2)
             W2 = 1.-W1
             K2 = BGRIDP(6*IBGRID)
             CALL SINTRP (W1, W2, BSPECS(1,1,K1,1), BSPECS(1,1,K2,1),
     &                    AC2(1,1,INDXGR), SPCDIR, SPCSIG)
!            --- store Hs from boundary condition in array HSOBND
             ETOT = 0.
             DO IS = 1, MSC
               SIG2 = SPCSIG(IS) ** 2
               DO ID = 1, MDC
                 ETOT = ETOT + SIG2 * AC2(ID,IS,INDXGR)
               ENDDO
             ENDDO
             IF (ETOT.GT.0.) THEN
               HS = 4. * SQRT(FRINTF*DDIR*ETOT)
             ELSE
               HS = 0.
             ENDIF
             COMPDA(INDXGR,JHSIBC) = HS
!
!            --- test output: parameters in test points on boundary
!
             IF (NPTST.GT.0) THEN
               DO IPTST = 1, NPTST
                 IF (OPTG.NE.5) THEN                                      40.80
                    IXP = XYTST(2*IPTST-1)
                    IYP = XYTST(2*IPTST)
                    LPB = INDXGR.EQ.KGRPNT(IXP,IYP)                       40.80
                 ELSE                                                     40.80
                    IVP = XYTST(IPTST)                                    40.80
                    LPB = INDXGR.EQ.IVP                                   40.80
                 ENDIF                                                    40.80
                 IF ( LPB ) THEN                                          40.80
                   IF (OPTG.NE.5) THEN                                    40.80
                      WRITE(PRTEST,72) IBGRID, IXP-1, IYP-1,
     &                                 W1, K1, W2, K2
  72                  FORMAT (' boundary point', 3I8, 2(F8.3, I4))
                   ELSE                                                   40.80
                      WRITE(PRTEST,73) IBGRID, IVP,                       40.80
     &                                 W1, K1, W2, K2                     40.80
  73                  FORMAT (' boundary vertex', 2I8, 2(F8.3, I4))       40.80
                   ENDIF                                                  40.80
                   AX = 0.
                   AY = 0.
                   ATOT = 0.
                   ASTOT = 0.
                   DO ID = 1, MDC
                     AADD = 0.
                     ASADD = 0.
                     DO IS = 1, MSC
                       SIG = SPCSIG(IS)
                       AA  = SIG*AC2(ID,IS,INDXGR)
                       AADD = AADD + AA
                       ASADD = ASADD + SIG*AA
                     ENDDO
                     AX = AX + AADD * SPCDIR(ID,2)
                     AY = AY + AADD * SPCDIR(ID,3)
                     ATOT = ATOT + AADD
                     ASTOT = ASTOT + ASADD
                   ENDDO
                   IF (ASTOT.GT.0.) THEN
                     HS = 4. * SQRT(FRINTF*DDIR*ASTOT)
                     APER = PI2 * ATOT / ASTOT
                     ADEG = 180./PI * ATAN2(AY,AX)
                   ELSE
                     HS = 0.
                     APER = -999.
                     ADEG = -999.
                   ENDIF
                   WRITE (PRTEST, 74) HS, APER, ADEG
  74               FORMAT (' Hs, Per, Dir: ', 3E12.4)
                 END IF
               END DO
             END IF
           END IF
         END DO
      END IF
!
!     --- update input fields (wind, water level etc.)
!
!     fields 5 and 6: wind
      IF (IFLDYN(5) .EQ. 1) THEN
        CALL FLFILE ( 5, 6, WXI, WYI,                                     40.31
     &               0, JWX2, JWX3, 0, JWY2, JWY3,
     &               COSWC, SINWC,                                        40.31 30.90
     &               COMPDA, XCGRID, YCGRID,
     &               KGRPNT, IERR)
        IF (STPNOW()) RETURN                                              34.01
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' wind from file in test points'
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JWX3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JWX3)                          40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' X-comp: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JWY3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JWY3)                          40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' Y-comp: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
 122      FORMAT (A, 10(1X,E11.4))
        ENDIF
      ENDIF

!     field 4: friction coeff.
      IF (IFLDYN(4) .EQ. 1) THEN
        CALL FLFILE ( 4, 0, FRIC, 0.,                                     40.31
     &               0, JFRC2, JFRC3, 0, 0, 0,
     &               1., 0.,                                              40.31 30.90
     &               COMPDA, XCGRID, YCGRID,
     &               KGRPNT, IERR)
        IF (STPNOW()) RETURN                                              34.01
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' fric coeff from file in test points'
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JFRC3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JFRC3)                         40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) 'friction: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
      ENDIF

!     field 7: water level
      IF (IFLDYN(7) .EQ. 1) THEN
        CALL FLFILE ( 7, 0, WLEVL, 0.,                                    40.31
     &               JWLV1, JWLV2, JWLV3, 0, 0, 0,
     &               1., 0.,                                              40.31 30.90
     &               COMPDA, XCGRID, YCGRID,
     &               KGRPNT, IERR)
        IF (STPNOW()) RETURN                                              34.01
!       Add bottom level to obtain depth
!       structured grid                                                   40.80
        DO 31 IX = 1, MXC
          DO 41 IY = 1, MYC
            INDX = KGRPNT(IX,IY)
            IF (INDX.GT.1) THEN
              XP = XCGRID(IX,IY)
              YP = YCGRID(IX,IY)
              DEP = SVALQI (XP, YP, 1, DEPTH, 1 ,IX ,IY)                  40.31 30.90
              COMPDA(INDX,JDP1) = COMPDA(INDX,JDP2)
              WLVL = COMPDA(INDX,JWLV2)
              DEPW = DEP + WLVL + WLEV
              COMPDA(INDX,JDP2) = DEPW
              IF (LSETUP.GT.0) THEN                                       40.14
                COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)                   40.14
              ENDIF                                                       40.14
            ENDIF
 41       CONTINUE
 31     CONTINUE
!       unstructured grid                                                 40.80
        DO INDX = 1, nverts                                               40.80
           XP = xcugrd(INDX)
           YP = ycugrd(INDX)
           IF ( IGTYPE(1).EQ.3 ) THEN
              DEP = DEPTH(INDX)
           ELSE
              DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)
           ENDIF
           COMPDA(INDX,JDP1) = COMPDA(INDX,JDP2)
           WLVL = COMPDA(INDX,JWLV2)
           DEPW = DEP + WLVL + WLEV
           COMPDA(INDX,JDP2) = DEPW
           IF (LSETUP.GT.0) THEN
              COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)
           ENDIF
        ENDDO                                                             40.80
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' water level from file in test points'
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JWLV3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JWLV3)                         40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' W-level: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
      ENDIF

!     field 2 and 3: current velocity
      IF (IFLDYN(2) .EQ. 1) THEN
        CALL FLFILE ( 2, 3, UXB, UYB,                                     40.31
     &               JVX1, JVX2, JVX3, JVY1, JVY2, JVY3,
     &               COSVC, SINVC,                                        40.31 30.90
     &               COMPDA, XCGRID, YCGRID,
     &               KGRPNT, IERR)
        IF (STPNOW()) RETURN                                              34.01
!       reduce current velocity if Froude number is larger than PNUMS(18)
!       structured grid                                                   40.80
        DO IX = 1, MXC                                                    40.13
          DO IY = 1, MYC                                                  40.13
            INDX = KGRPNT(IX,IY)                                          40.13
            IF (INDX.GT.1) THEN                                           40.13
              DEPW = COMPDA(INDX,JDP2)
              IF (DEPW.GT.0.) THEN
                UU = COMPDA(INDX,JVX2)
                VV = COMPDA(INDX,JVY2)
                VTOT = SQRT (UU*UU + VV*VV)
                CGMAX = PNUMS(18)*SQRT(GRAV*DEPW)
                IF (VTOT .GT. CGMAX) THEN
                  CGFACT = CGMAX / VTOT
                  COMPDA(INDX,JVX2) = UU * CGFACT
                  COMPDA(INDX,JVY2) = VV * CGFACT
!                 write IX,IY to error points file
                  IF (ERRPTS.GT.0.AND.IAMMASTER) THEN                     40.95 40.30
                    WRITE (ERRPTS, 211) IX+MXF-1, IY+MYF-1, 1
 211                FORMAT (I4, 1X, I4, 1X, I2)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF                                                         40.13
          ENDDO                                                           40.13
        ENDDO                                                             40.13
!       unstructured grid                                                 40.80
        DO INDX = 1, nverts                                               40.80
           DEPW = COMPDA(INDX,JDP2)
           IF (DEPW.GT.0.) THEN
              UU = COMPDA(INDX,JVX2)
              VV = COMPDA(INDX,JVY2)
              VTOT = SQRT (UU*UU + VV*VV)
              CGMAX = PNUMS(18)*SQRT(GRAV*DEPW)
              IF (VTOT .GT. CGMAX) THEN
                 CGFACT = CGMAX / VTOT
                 COMPDA(INDX,JVX2) = UU * CGFACT
                 COMPDA(INDX,JVY2) = VV * CGFACT
!                write INDX to error points file
                 IF (ERRPTS.GT.0) THEN
                    WRITE (ERRPTS, 212) INDX, 1
 212                FORMAT (I4, 1X, I2)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO                                                             40.80
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' current vel from file in test points'
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JVX3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JVX3)                          40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' X-comp: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
                IXP = XYTST(2*IPTST-1)
                IYP = XYTST(2*IPTST)
                TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JVY3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JVY3)                          40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' Y-comp: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
      ENDIF

!     field 11: field containing number of plants per square meter
      IF (IFLDYN(11) .EQ. 1) THEN
        CALL FLFILE (11, 0, NPLAF, 0.,
     &               0, JNPLA2, JNPLA3, 0, 0, 0,
     &               1., 0.,
     &               COMPDA, XCGRID, YCGRID,
     &               KGRPNT, IERR)
        IF (STPNOW()) RETURN
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' # plants/m2 from file in test points'
          IF (OPTG.NE.5) THEN                                             40.80
             DO IPTST = 1, MIN(10,NPTST)
               IXP = XYTST(2*IPTST-1)
               IYP = XYTST(2*IPTST)
               TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JNPLA3)
             ENDDO
          ELSE                                                            40.80
             DO IPTST = 1, MIN(10,NPTST)                                  40.80
                IVP = XYTST(IPTST)                                        40.80
                TSTVAL(IPTST) = COMPDA(IVP,JNPLA3)                        40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
          WRITE (PRTEST, 122) ' veg dens: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
      ENDIF
!
!     End of subroutine SNEXTI
!
      RETURN
      END
!****************************************************************
!
      SUBROUTINE RBFILE (SPCSIG, SPCDIR, BFILED, BSPLOC,
     &                   BSPDIR, BSPFRQ, BSPECS, XYTST )                  40.31 30.90
!
!****************************************************************
      USE TIMECOMM                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      IMPLICIT NONE

      REAL,   INTENT(IN)     ::  SPCDIR(MDC,6)                            30.82
      REAL,   INTENT(IN)     ::  SPCSIG(MSC)                              30.82
      REAL,   INTENT(INOUT)  ::  BSPECS(MDC,MSC,NBSPEC,2)
      REAL   , INTENT(INOUT)  :: BSPDIR(*)
      REAL   , INTENT(INOUT)  :: BSPFRQ(*)
      INTEGER, INTENT(INOUT)  ::  BFILED(*)
      INTEGER, INTENT(INOUT)  ::  BSPLOC(*)
      INTEGER, INTENT(IN)     ::  XYTST(*)
      INTEGER   DORDER, NDSL, NDSD, IHD, IBOUNC, IBSPEC, IERR
      INTEGER   ID, IS, JJ, NANG, NFRE, COUNT_IT, IENT, II               40.05
      INTEGER   WWDATE, WWTIME                                           40.05
      REAL      TIMF1, TIMF2, UFAC, W1, BDEPTH, DUM_A, RFAC               40.05
      REAL      XLON, XLAT, XDATE, EMEAN, THQ, FMEAN, USNEW, THWNEW

      DOUBLE PRECISION DDATE
!
      LOGICAL   NSTATF, UNFORM
!
      CHARACTER BTYPE *4, HEDLIN *80, TIMSTR *20,
     &          DATITM(5) *18, PTNME*12
      CHARACTER (LEN=14) :: CDATE     ! date-time
!
      REAL, ALLOCATABLE :: SPAUX(:,:)

      SAVE COUNT_IT                                                       40.31
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'RBFILE')
!
!
!     if data file is exhausted, return
      IF (BFILED(1).EQ.-1) GOTO 900
      NSTATF = (BFILED(1) .GT. 0)
      IERR   = 0
      CALL COPYCH (BTYPE, 'F', BFILED(7), 1, IERR)
!
      IF (BTYPE.EQ.'WAMW' .OR. BTYPE.EQ.'WAMC') THEN
        UNFORM = .TRUE.
      ELSE
        UNFORM = .FALSE.
      ENDIF
!
      IF (BFILED(2).LT.0) COUNT_IT = 0                                    40.05
!
      DORDER = BFILED(9)
      TIMF1 = REAL(BFILED(2))
      TIMF2 = REAL(BFILED(3))
      IF (ITEST.GE.120) WRITE (PRINTF, 187) BFILED(1), BTYPE,
     &      TIMF1, TIMF2, TIMCO
 187  FORMAT (' Boundary', I2, 2X, A, ' times: ', 3F10.1)
!
      ALLOCATE(SPAUX(MDC,MSC))                                            40.31
!
!     if present time > time of last set of spectra, read new spectra
!
!     While Loop named GLOOP
!
      GLOOP : DO                                                          40.05
!
        IF (TIMCO.GT.TIMF2) THEN                                          40.05
!
!     then read from boundary nesting files all the information
!     and the spectral ones
!
!     COUNT_IT : counter of the time which enter in a WW3  boundary file
!     and read spectral - used to calculate NHED (number of heading lines
!     per time) in WW3N case.
!
          COUNT_IT = COUNT_IT+1                                           40.05
          NDSL = BFILED(4)
          NDSD = BFILED(5)
!
!     in WW3 case rewind the boundary file
          IF(BTYPE.EQ.'WW3N') REWIND(NDSD)                                40.05
!
          TIMF1 = TIMF2
!
!         move new spectra to old for all boundary points
          DO IBOUNC = 1, BFILED(8)
            IBSPEC = BSPLOC(IBOUNC)

            DO ID = 1, MDC
              DO IS = 1, MSC
                BSPECS(ID,IS,IBSPEC,1) = BSPECS(ID,IS,IBSPEC,2)
              ENDDO
            ENDDO

            IF (ITEST.GE.80) WRITE (PRTEST, *) ' spectrum moved ',
     &         IBSPEC, TIMF1
          ENDDO
!
 210       CONTINUE
!
!     define the new  NHED in every time step for WW3 case
          IF (BTYPE.EQ.'WW3N') THEN
            BFILED(15) = BFILED(14)+(1+((BFILED(16)-1)+
     &                 CEILING((BFILED(12)* BFILED(10))/7.))*
     &                 BFILED(8))*(COUNT_IT-1)                            40.05
          ENDIF                                                           40.05
!
!     read heading lines per time step
!
!     HBFL is loop over the number of heading lines per time step
          HBFL : DO IHD = 1, BFILED(15)                                   40.05
            IF (UNFORM) THEN
              READ (NDSD, END=380, ERR=920)
            ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN             40.05
              READ (NDSD,*)                                               40.05
            ELSE                                                          40.05
              READ (NDSD, '(A)', END=380, ERR=920) HEDLIN
              IF (ITEST.GE.90) WRITE (PRINTF, 212) HEDLIN
 212          FORMAT (' heading line: ', A)
              IF (BTYPE.EQ.'SWNT') THEN
!               convert time string to time in seconds
                CALL DTRETI (HEDLIN(1:18), BFILED(6), TIMF2)              40.03
              ENDIF
            ENDIF
          ENDDO HBFL
!
          IF (.NOT.NSTATF) TIMF2 = 0.
          IF (ITEST.GE.60) WRITE (PRINTF, 214)
     &        TIMF1, TIMF2, TIMCO, BFILED(8), BFILED(15)
 214      FORMAT (' Boundary times ', 3F12.0, 2X, 4I4)
!
!         read additional information from the headers and spectrum from
!         the nesting boundary files (for all the cases)
!
!
!       BP_LOOP loop over the boundary nesting points
          BP_LOOP : DO  IBOUNC = 1, BFILED(8)                             40.05
!
            IBSPEC = BSPLOC(IBOUNC)
!           division by 2*PI to account for difference in definition of freq.
!           Hz to rad/s
            NANG  = BFILED(10)                                            40.00
!
!           calculate UFAC for the different nesting cases
            IF (BTYPE(1:3).EQ.'SWN' .AND. NANG.GT.0) THEN
!           in addition multiply by 180/PI to account for directions in degr
!           instead of radians (SWAN 2D spectral files)
              UFAC = 180./ (2.*PI**2)
            ELSE                                                          40.05
              UFAC = 1./ (2.*PI)
            ENDIF
!           divide by Rho*Grav if quantity in file is energy density
            IF ((BFILED(17).EQ.1))  UFAC = UFAC / (RHO*GRAV)
!
!           read information from heading lines per spectrum
!
!          do loop over the numbers of the heading lines per spectrum
!
            DO IHD = 1, BFILED(16)
!
              IF (IBOUNC.EQ.1) THEN                                       40.03
!             for the first spectrum (first point), read time from
!             heading line                                                40.03
                IF (BTYPE.EQ.'WAMW') THEN                                 40.03
                  READ(NDSD, END=380, ERR=920) XLON, XLAT,
     &               CDATE(1:LWDATE),                                     41.13
     &               EMEAN, THQ, FMEAN
                  IF (LEN_TRIM(CDATE) == 10 ) THEN
                     TIMSTR = TRIM(CDATE)
                  ELSEIF (LEN_TRIM(CDATE) == 12 ) THEN
                     TIMSTR = CDATE(1:10)
                  ELSEIF (LEN_TRIM(CDATE) == 14 ) THEN
                     TIMSTR = CDATE(3:12)
                  ENDIF
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  40.03 40.05
                ELSE IF (BTYPE.EQ.'WAMC') THEN
                  READ(NDSD, END=380, ERR=920) XLON, XLAT, XDATE,
     &                EMEAN, THQ, FMEAN, USNEW, THWNEW
                  WRITE (TIMSTR,'(F11.0,9X)') XDATE
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  40.03 40.05
                ELSE IF (BTYPE.EQ.'WAMF') THEN
                  READ(NDSD,*, END=380, ERR=920) XLON, XLAT, DDATE,
     &              EMEAN, THQ, FMEAN, USNEW, THWNEW
                  WRITE (TIMSTR,'(F11.0,9X)') DDATE
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  40.03 40.05
                ELSE IF (BTYPE.EQ.'WW3N') THEN                            40.05
!                  read from heading lines per spectrum , date ,time      40.05
                  IF(IHD.EQ.1) THEN                                       40.05
                    READ(NDSD,*, END=380, ERR=920) WWDATE, WWTIME         40.05
                    WRITE (TIMSTR, 118) WWDATE, WWTIME                    40.05
 118                FORMAT (I8,'.',I6)                                    40.05
!                    convert time string to time in seconds
                    CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                40.05
                  ELSE                                                    40.05
!                   DUM_A dummy local variable used to read formatted     40.05
!                   files                                                 40.05
!                   this variable will be not used in any calculation     40.05
                    READ (NDSD,901) PTNME, DUM_A, DUM_A, BDEPTH,          40.05
     &                              DUM_A, DUM_A, DUM_A, DUM_A            40.05
                  ENDIF                                                   40.05
                ELSE                                                      40.05
!                 SWAN files                                              40.05
                  READ (NDSD, '(A)', END=380, ERR=920) HEDLIN             40.05
                  IF (ITEST.GE.100) WRITE (PRINTF, 212) HEDLIN            40.05
                ENDIF
              ELSE
                IF (UNFORM) THEN
                  READ (NDSD, END=380, ERR=920)
                ELSE IF (BTYPE.EQ.'WW3N') THEN                            40.13 40.05
!                 read from heading lines per spectrum depth of the       40.05
!                 boundary point                                          40.05
                  IF (IHD.GT.1) EXIT                                      40.05
!                 exit is used because the header per spectrum in WW3     40.05
!                 for IBOUNC>1 is one line                                40.05
!                 DUM_A is a dummy local parameter used in formatted read 40.05
                  READ (NDSD,901) PTNME, DUM_A, DUM_A, BDEPTH,
     &                            DUM_A, DUM_A, DUM_A, DUM_A              40.05
                ELSE IF (BTYPE.EQ.'WAMF') THEN
!                 read HEDLIN replaced because data are sometimes written 40.13
!                 on two subsequent lines                                 40.13
                  READ(NDSD,*, END=380, ERR=920) XLON, XLAT, DDATE,       40.13
     &              EMEAN, THQ, FMEAN, USNEW, THWNEW
                ELSE
                  READ (NDSD, '(A)', END=380, ERR=920) HEDLIN
                  IF (ITEST.GE.100) WRITE (PRINTF, 212) HEDLIN
                ENDIF
              ENDIF
!
              IF (BTYPE(1:3).EQ.'SWN') THEN
!             SWAN nesting: take proper action if heading line contains ZERO or NODATA
                IF (HEDLIN(1:6).EQ.'NODATA' .OR. HEDLIN(1:4).EQ.'ZERO')
     &               THEN
                  DO IS = 1, MSC
                    DO ID = 1, MDC
                      BSPECS(ID,IS,IBSPEC,2) = 0.
                    ENDDO
                  ENDDO
!                 skip reading of values
                  CYCLE BP_LOOP                                           40.05
                ELSE IF (HEDLIN(1:6).EQ.'FACTOR') THEN
                  READ (NDSD, *) RFAC
                  UFAC = UFAC * RFAC
!                 multiply factor read from file by UFAC (factor following from
!                 type of file)
                ELSE
!                 note: in case of 1D spectra heading line can be ignored
                  IF (NANG.GT.0) THEN                                     40.00
                    CALL MSGERR (3,
     &                'incorrect code in b.c. file: '//HEDLIN(1:20))      40.00
                  ENDIF                                                   40.00
                ENDIF
              ENDIF
!           end loop over heading lines per spectrum
            ENDDO
!
!           test output: which spectrum is processed                      40.00
!
            IF (ITEST.GE.60) THEN
              INQUIRE (UNIT=NDSD, NAME=FILENM)                            40.00
              WRITE (PRTEST, 188) FILENM, CHTIME, IBOUNC, UFAC
 188          FORMAT
     &   (' read spectrum ', A, '; time=', A, ' nr=', I3, F9.3)           40.00
            ENDIF
!
!       start reading incoming wave data
!
            IF (BTYPE.EQ.'TPAR') THEN
              READ (NDSD, 222, END=380) HEDLIN
 222          FORMAT (A)
              CALL LSPLIT (HEDLIN, DATITM, 5)
              CALL DTRETI (DATITM(1), BFILED(6), TIMF2)                   40.03
              DO II = 1, 4
                READ (DATITM(II+1), '(G12.0)') SPPARM(II)
              ENDDO
              IF (ITEST.GE.60) WRITE (PRTEST, *) ' TPAR boundary ',
     &              TIMF2, (SPPARM(JJ), JJ=1,4)
              BFILED(3) = NINT(TIMF2)
              CALL SSHAPE (BSPECS(1,1,IBSPEC,2), SPCSIG, SPCDIR,
     &                 FSHAPE, DSHAPE)
            ELSE
!         other (spectral) boundary conditions
              NANG  = BFILED(10)
              NFRE  = BFILED(12)
!
!         call RESPEC subroutine to read the spectum of the bound. files
!
             CALL RESPEC (BTYPE, NDSD, BFILED, UNFORM, DORDER,            40.00
     &           SPCSIG, SPCDIR, BSPFRQ, BSPDIR, BSPECS(1,1,IBSPEC,2),    40.31 30.90
     &           UFAC, IERR)
             IF (IERR.EQ.9) GOTO 380
            ENDIF
!
          END DO BP_LOOP                                                  40.05
!
          IF (ITEST.GE.60) WRITE (PRINTF, 287) BTYPE, TIMF2
 287      FORMAT
     &         (' Boundary data type ', A, ' processed, time: ', F10.1)
!
!        cycle back to will loop
          CYCLE GLOOP                                                     40.05
!
!         if there are no more data on a boundary data file
!         close this file, and see if there is a next one
!
 380      CLOSE(NDSD)
!         read filename of next boundary file and open them
          IF (NDSL.GT.0) THEN                                             40.05
            READ (NDSL, '(A)', END=390, ERR=930) FILENM
            IF (UNFORM) THEN
              OPEN (NDSD, FILE=FILENM, FORM='UNFORMATTED',
     &               STATUS='OLD', ERR=930)
!              read heading lines
              DO IHD = 1, BFILED(14)                                      40.03
                READ (NDSD, END=940, ERR=920)
              ENDDO
            ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN             40.05
!             if it is WW3 open the new file only
              OPEN (NDSD, FILE=FILENM, FORM='FORMATTED',                  40.05
     &               STATUS='OLD', ERR=930)                               40.05
              COUNT_IT = 1                                                40.05
            ELSE                                                          40.05
              OPEN (NDSD, FILE=FILENM, FORM='FORMATTED',
     &               STATUS='OLD', ERR=930)
              DO IHD = 1, BFILED(14)                                      40.03
!               read heading lines
                READ (NDSD, '(A)', END=940, ERR=920) HEDLIN
                IF (ITEST.GE.80) WRITE (PRINTF, 212) HEDLIN
              ENDDO
            ENDIF
!
!      go back to statement 210 to start again the procedure of reading
!      info and spectrum from the new boundary file
!
            GOTO 210
!
          ELSE
!           boundary data are read from a single file
            GOTO 392                                                      40.13
          ENDIF                                                           40.05
!         close file containing filenames
 390      CLOSE (NDSL)
!
          BFILED(5) = 0
!         write message and close file containing spectra
 392      CALL MSGERR (1, 'data on boundary file exhausted')
!
          BFILED(4) = 0
          BFILED(1) = -1
          TIMF2 = 999999999.
          CYCLE GLOOP                                                     40.05
!
!        (if necessary) data have been read from file, now interpolate in time
!
        ELSE                                                              40.05
!       if present time <= time of last set of spectra then
!       transform to spectral resolution used in SWAN to
!       obtain boundary spectra
!
          IF (TIMF1.NE.TIMF2) THEN                                        40.41
             W1 = (TIMF2-TIMCO) / (TIMF2-TIMF1)
          ELSE
             W1 = 0.
          END IF
          DO IBOUNC = 1, BFILED(8)
            IBSPEC = BSPLOC(IBOUNC)
            IF (IBOUNC.EQ.1 .AND. ITEST.GE.80) WRITE (PRTEST, 403)
     &           TIMCO, W1, TIMF1, TIMF2, IBSPEC
 403        FORMAT (' interp in time ', F14.1, F8.3, 2F14.1, I4)
!
!       interpolate spectra in time; result has to be store in BSPECS(..,1)
!       first interpolate to auxiliary array
!
            CALL SINTRP (W1, 1.-W1, BSPECS(1,1,IBSPEC,1),
     &               BSPECS(1,1,IBSPEC,2), SPAUX,                         40.31 30.90
     &               SPCDIR, SPCSIG)
!       use SINTRP to copy contents of aux. array to BSPECS(..,1)
!
            CALL SINTRP (1., 0., SPAUX,                                   40.31 30.90
     &               BSPECS(1,1,IBSPEC,2), BSPECS(1,1,IBSPEC,1),
     &               SPCDIR, SPCSIG)
          ENDDO
          BFILED(2) = NINT(TIMCO)
          BFILED(3) = NINT(TIMF2)
          EXIT GLOOP
!       end of time comparison
        ENDIF                                                             40.05

!     end of while loop
      END DO GLOOP                                                        40.05
!
!
 901  FORMAT (A12,2F7.2,F10.1,2(F7.2,F6.1))
!
      DEALLOCATE(SPAUX)                                                   40.31
!
 900  RETURN
!
 920  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    40.00
      CALL MSGERR (4,
     &     'error reading data from boundary file '//FILENM)              40.00
      GOTO 900
 930  CALL MSGERR (4,
     &     'error opening boundary file '//FILENM)                        40.00
      GOTO 900
 940  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    40.00
      CALL MSGERR (4,
     &     'unexpected end of file on boundary file '//FILENM)            40.00
      GOTO 900
!
!     End of subroutine RBFILE
!
      END
!****************************************************************
!
      SUBROUTINE RESPEC (BTYPE, NDSD, BFILED, UNFORM, DORDER,             40.00
     &                   SPCSIG, SPCDIR, BSPFRQ, BSPDIR, LSPEC, UFAC,
     &                   IERR)
!
!****************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
!
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)  :: BFILED(*)                                40.05
!
      REAL,    INTENT(IN)     :: SPCDIR(MDC,6)                            30.82
      REAL,    INTENT(IN)     :: SPCSIG(MSC)                              30.82
      REAL,    INTENT(IN)     :: BSPDIR(*)
      REAL,    INTENT(IN)     :: BSPFRQ(*)
      REAL,    INTENT(INOUT)  :: LSPEC(MDC,MSC)
      REAL,    INTENT(IN)     :: UFAC                                     40.05
!
      INTEGER,  INTENT(IN) ::  NDSD                                       40.00
      INTEGER,  INTENT(IN) ::  DORDER                                     40.00
      INTEGER,  INTENT(INOUT) ::  IERR                                    40.00
!
      LOGICAL   UNFORM
!
      CHARACTER BTYPE *4
!
      INTEGER   IANG, IFRE, ID, IS, ISIGTA,IENT,NFRE,NANG
!
      REAL, ALLOCATABLE :: BAUX0(:,:)                                     40.05
      REAL, ALLOCATABLE :: BAUX1(:,:)                                     40.31
      REAL, ALLOCATABLE :: BAUX2(:,:)                                     40.31
      REAL, ALLOCATABLE :: BAUX3(:)                                       40.31
      REAL, ALLOCATABLE :: BAUX4(:)                                       40.31

      REAL      ETOT, ADEG, DD, ADIR, MS, CTOT, ACOS, CDIR
      REAL      DSUM, DSPR

      REAL :: GAMMA                                                       40.03
      LOGICAL EQREAL                                                      40.41

      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'RESPEC')
!
!     read the values from array BFILED for the number of
!     direction and frequencies
!
      NANG = BFILED(10)                                                   40.05
      NFRE = BFILED(12)                                                   40.05
      ALLOCATE (BAUX0(NFRE,NANG))                                         40.02
      ALLOCATE (BAUX1(NANG,NFRE))                                         40.31
      ALLOCATE (BAUX2(MDC ,NFRE))                                         40.31
      ALLOCATE (BAUX3(NFRE))                                              40.31
      ALLOCATE (BAUX4(MSC))                                               40.31

      IF (ITEST.GE.60) THEN
        INQUIRE (UNIT=NDSD, NAME=FILENM)                                  40.00
        WRITE (PRTEST, 7) FILENM, NANG, NFRE, UFAC, BFILED(18),           40.00
     &  BFILED(19), UNFORM
   7    FORMAT (' Entry RESPEC, reading file:', A, /, 6X,
     &  I3, ' angles ', I3, ' freqs; factor ',
     &  E12.4, ' dir.def:', I2, ' dir.spr.def:', I2, ' unform:', L2)      40.00
      ENDIF
!
!     Initialisation
!
      ISIGTA = MSC                                                        40.02
!
!     read spectral energy densities from b.c. file
      IF (NANG.EQ.0) THEN
!
!     1-D spectral input                                                  40.05
!
        DO IFRE=1,NFRE                                                    40.05
!         1-D spectral input (only function of frequency)
          IF (IFRE.EQ.1) THEN
            READ (NDSD,*,END=940,ERR=920) ETOT, ADEG, DD
          ELSE
            READ (NDSD,*,END=930,ERR=920) ETOT, ADEG, DD
          ENDIF
!
          IF (EQREAL(ETOT,REAL(BFILED(11)))) THEN                         40.41
!            in case of exception value, no energy
             ETOT = 0.
             ADEG = 0.
             DD   = 180./PI
          END IF
!
          IF (BFILED(18).EQ.1) THEN                                       40.00
!           conversion from degrees to radians
            ADIR = PI * ADEG / 180.
          ELSE
!           conversion from Nautical to Cartesian conv.
            ADIR = PI * (180.+DNORTH-ADEG) / 180.
          ENDIF
!
          IF (BFILED(19).EQ.1) THEN                                       40.00
!           DSPR is directional spread in radians
            DSPR = PI * DD / 180.
            IF (DSPR.NE.0.) THEN                                          40.41
               MS = MAX (DSPR**(-2) - 2., 1.)
            ELSE
               MS = 1000.                                                 40.41
            END IF
          ELSE
            MS = DD
          ENDIF
!
          IF (ITEST.GE.80) WRITE (PRTEST, 12) IFRE, ETOT,                 40.00
     &          180.*ADIR/PI, MS
  12      FORMAT (' read freq ', I3, E10.3, '; Cart dir ',                40.00
     &          F7.1, '; Cos power ', F7.2)
!
!         generate distribution over directions
!
!         equations taken from Jahnke & Emde (chapter Factorial Function)
          IF (MS.GT.10.) THEN
            CTOT = SQRT(MS/(2.*PI)) * (1. + 0.25/MS)                      30.81 40.00
          ELSE
            CTOT = 2.**MS * (GAMMA(1.+0.5*MS))**2 / (PI * GAMMA(1.+MS))   40.00
          ENDIF
          DSUM = 0.
          DO ID = 1, MDC
            ACOS = COS(SPCDIR(ID,1) - ADIR)                               30.82
            IF (ACOS .GT. 0.) THEN
              CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
            ELSE
              CDIR = 1.E-10
            ENDIF
            IF (ITEST.GE.20) DSUM = DSUM + CDIR * DDIR                    40.00
            BAUX2(ID,IFRE) = CDIR * ETOT
          ENDDO
          IF (ITEST.GE.20) THEN
            IF (ABS(DSUM-1.).GT.0.1) WRITE (PRTEST, 138) DSUM, CTOT, MS   40.00
 138        FORMAT (' integral over directions is ', F9.4,
     &              ' with CTOT=', F10.3,'; power=', F8.2)                40.00
          ENDIF
        END DO                                                            40.05
!
      ELSE
!
!   (2D) fully spectral input
        IF (UNFORM) THEN
!          unformatted reading
           IF (DORDER.LT.0) THEN
              READ (NDSD,END=930,ERR=920)
     &               ((BAUX1(IANG,IFRE),IANG=NANG,1,-1),IFRE=1,NFRE)
           ELSE
              READ (NDSD,END=930,ERR=920)
     &               ((BAUX1(IANG,IFRE),IANG=1,NANG),IFRE=1,NFRE)
           ENDIF
        ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN                 40.05
!
!         WW3 reading
!         BAUX0 local array to read the energy spectra from boundary files
!
          READ (NDSD,902,END=940,ERR=920)                                 40.15
     &      ((BAUX0(IFRE,IANG),IFRE=1,NFRE),IANG=1,NANG,1)                40.15
 902      FORMAT (7E11.3)                                                 40.15
!
!         energy (variance) density   from E(FRQ,TH) to  E(TH,FRQ)
!
          DO IANG = 1,NANG                                                40.05
            DO IFRE = 1,NFRE                                              40.05
              BAUX1(IANG,IFRE) = BAUX0(IFRE,IANG)                         40.05
            ENDDO                                                         40.05
          ENDDO                                                           40.05
        ELSE
!         format reading (except WW3)
          IF (DORDER.LT.0) THEN                                           40.31
            READ (NDSD,*,END=930,ERR=920)                                 40.61
     &             ((BAUX1(IANG,IFRE),IANG=NANG,1,-1),IFRE=1,NFRE)        40.61
          ELSE                                                            40.31
            READ (NDSD,*,END=930,ERR=920)                                 40.61
     &             ((BAUX1(IANG,IFRE),IANG=1,NANG),IFRE=1,NFRE)           40.61
          ENDIF                                                           40.31
        ENDIF
!
        IF (ITEST.GE.120) THEN
          WRITE (PRINTF,*)' Spectra from file'
          DO IFRE = 1, NFRE
            WRITE (PRINTF,*) IFRE, (BAUX1(IANG,IFRE),IANG=1,NANG)
          ENDDO
        ENDIF

!       --- in case of exception value, no energy                         40.41

        DO IANG = 1,NANG
           DO IFRE = 1,NFRE
              IF (EQREAL(BAUX1(IANG,IFRE),REAL(BFILED(11)))) THEN
                 BAUX1(IANG,IFRE) = 0.
              END IF
           END DO
        END DO
!
!       --- transform to spectral directions used in SWAN
!           results appear in array BAUX2(MDC,NFRE)
!
        DO IFRE = 1, NFRE                                                 40.05
          CALL CHGBAS (BSPDIR, SPCDIR, PI2, BAUX1(1,IFRE),
     &                 BAUX2(1,IFRE), NANG, MDC, ITEST, PRTEST)
        END DO                                                            40.05
      ENDIF
!
!     interpolate energy densities to SWAN frequencies distribution
!
      IF (BSPFRQ(NFRE) .LT. SPCSIG(MSC)) THEN
        DO IS = MSC, 1, -1
          IF (SPCSIG(IS).LT.BSPFRQ(NFRE)) THEN
!           ISIGTA is the last frequency which is determined by interpolation
!           higher frequencies are determined by tail expression
            ISIGTA = IS
            EXIT                                                          40.05
          ENDIF
        ENDDO
      ELSE
        ISIGTA = MSC
      ENDIF
!
!
!     UFAC is the product of the multiplication factor read from file
!     and the factor to transform from energy/Hz to energy/(rad/s)
!     and from energy/degr to energy/rad (latter only for 2d spectra)     40.00
!
      DO  ID = 1,MDC                                                      40.05
        DO  IFRE = 1,NFRE                                                 40.05
          BAUX3(IFRE) = UFAC * BAUX2(ID,IFRE)
        ENDDO                                                             40.05

!       interpolate over frequency keeping energy constant, output BAUX4(MSC)
        CALL CHGBAS (BSPFRQ, SPCSIG, 0., BAUX3, BAUX4, NFRE, MSC,
     &               ITEST, PRTEST)
        DO  IS=1,MSC                                                      40.05
          IF (IS.LE.ISIGTA) THEN
!
!           to convert energy density to action density
!
            LSPEC(ID,IS) = BAUX4(IS)/SPCSIG(IS)
          ELSE
!
!           add a tail when IS > ISIGTA
!
            LSPEC(ID,IS) = LSPEC(ID,ISIGTA) *
     &                 (SPCSIG(ISIGTA)/SPCSIG(IS))**(PWTAIL(1)+1)
          ENDIF
          IF (ITEST.GE.140) THEN
            WRITE (PRTEST, *) 'ID,IS,LSPEC(ID,IS)', ID,IS,LSPEC(ID,IS)
          ENDIF
        ENDDO                                                             40.05
      ENDDO                                                               40.05
!
      DEALLOCATE (BAUX0,BAUX1,BAUX2,BAUX3,BAUX4)                          40.31 40.02

 900  IERR = 0
      RETURN
 920  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    40.00
      CALL MSGERR (2,
     &      'read error in boundary condition file '//FILENM)             40.00


      RETURN
 930  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    40.00
      CALL MSGERR (2,
     &      'insufficient data in boundary condition file '//FILENM)      40.00

 940  IERR = 9

      RETURN

      END SUBROUTINE RESPEC
!**********************************************************************
!
      SUBROUTINE FLFILE (IGR1, IGR2,
     &                   ARR, ARR2, JX1, JX2, JX3, JY1, JY2, JY3,         40.31
     &                   COSFC, SINFC, COMPDA,                            40.31 30.90
     &                   XCGRID, YCGRID,
     &                   KGRPNT, IERR)
!
!**********************************************************************

      USE TIMECOMM                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80

      IMPLICIT NONE

      LOGICAL STPNOW                                                      34.01

      INTEGER    KGRPNT(MXC,MYC),
     &           IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
!
      REAL       COMPDA(MCGRD,MCMVAR),
     &           XCGRID(MXC,MYC), YCGRID(MXC,MYC),
     &           COSFC, SINFC
      REAL       ARR(*), ARR2(*)
!
!     local variables
!
      INTEGER    IENT, INDX, IX, IY
!     INDX       counter of comp. grid points
!     IX         index in x-dir of comput grid point
!     IY         index in y-dir of comput grid point
!
      REAL       SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
      REAL       TIMR1, XP, YP, UU, VV, VTOT, W1, W3,
     &           SIZE1, SIZE2, SIZE3
!     TIMR1      time of one but last input field
!     XP         x-coord of one comput grid point
!     YP         y-coord of one comput grid point
!     UU         x-component of vector, or scalar value
!     VV         y-component of vector
!     VTOT       length of vector
!     W1         weighting coeff for interpolation in time
!     W3         weighting coeff for interpolation in time
!     DIRE       direction of interpolated vector
!     SIZE1      length of vector at time TIMR1
!     SIZE2      length of vector at time TIMCO
!     SIZE3      length of vector at time TIMR2
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'FLFILE')
!
      IERR = 0
!
      IF (JX1.GT.1) THEN
        DO INDX = 1, MCGRD
          COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
        ENDDO
      ENDIF
      IF (IGR2.GT.0 .AND. JY1.GT.1) THEN
        DO INDX = 1, MCGRD
          COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
        ENDDO
      ENDIF
      TIMR1 = TIMCO - DT
!
 200  IF (TIMCO.LE.IFLTIM(IGR1)) GOTO 400
      TIMR1 = IFLTIM(IGR1)
      IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
      IF (IFLTIM(IGR1) .GT. IFLEND(IGR1)) THEN
        IFLTIM(IGR1) = 1.E10
        IF (IGR2.GT.0) IFLTIM(IGR2) = IFLTIM(IGR1)
        GOTO 400
      ENDIF
      IF (IFLNDS(IGR1).GT.0) THEN                                         40.03
        IF (INODE.EQ.MASTER) THEN
!ADC!          --- if we have coupled to the quantities from ADCIRC,          41.20
!ADC!              then grab them from memory instead of reading them
!ADC!              from the external file
!ADC           IF ( (IGR1.EQ.7).AND.COUPWLV ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_ETA2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_ETA2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.2).AND.COUPCUR ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_UU2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_UU2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.5).AND.COUPWIND ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_WX2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_WX2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.4).AND.COUPFRIC ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_Z0(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_Z0(INDX,2))
!ADC              ENDDO
!ADC           ELSE
           CALL INAR2D( ARR, MXG(IGR1), MYG(IGR1),                        40.31 40.02
     &                  IFLNDF(IGR1),
     &                  IFLNDS(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),
     &                  IFLIDL(IGR1), IFLFAC(IGR1),
     &                  IFLNHD(IGR1), IFLNHF(IGR1))
           IF (STPNOW()) RETURN                                           34.01
!ADC           ENDIF
        END IF
        CALL SWBROADC(IFLIDL(IGR1),1,SWINT)                               40.30
        IF (IFLIDL(IGR1).LT.0) THEN
!         end of file was encountered
          IFLTIM(IGR1) = 1.E10
          IF (IGR2.GT.0) IFLTIM(IGR2) = IFLTIM(IGR1)
          GOTO 400
        ELSE
          CALL SWBROADC(ARR,MXG(IGR1)*MYG(IGR1),SWREAL)                   40.31 40.30
        ENDIF
      ELSE                                                                40.13
        IF (ITEST.GE.20) THEN                                             40.13
          CALL MSGERR (1,
     &    'no read of input field because unit nr=0')                     40.13
          WRITE (PRINTF, 208) IGR1
 208      FORMAT (' field nr.', I2)
        ENDIF
      ENDIF                                                               40.03
      IF (IGR2.GT.0) THEN
        IFLTIM(IGR2) = IFLTIM(IGR1)
        IF (IFLNDS(IGR2).GT.0) THEN                                       40.03
          IF (INODE.EQ.MASTER) THEN                                       40.30
!ADC!            added these lines to grab the y-components from memory       41.20
!ADC             IF ( (IGR2.EQ.3).AND.COUPCUR ) THEN
!ADC               DO INDX = 1, nverts
!ADC                 ARR2(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                   * REAL(ADCIRC_VV2(INDX,1))
!ADC     &                   + REAL(InterpoWeight)
!ADC     &                   * REAL(ADCIRC_VV2(INDX,2))
!ADC               ENDDO
!ADC             ELSEIF ( (IGR2.EQ.6).AND.COUPWIND ) THEN
!ADC               DO INDX = 1, nverts
!ADC                 ARR2(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                   * REAL(ADCIRC_WY2(INDX,1))
!ADC     &                   + REAL(InterpoWeight)
!ADC     &                   * REAL(ADCIRC_WY2(INDX,2))
!ADC               ENDDO
!ADC             ELSE
             CALL INAR2D( ARR2, MXG(IGR2), MYG(IGR2),                     40.31 40.02
     &                    IFLNDF(IGR2),
     &                    IFLNDS(IGR2), IFLIFM(IGR2), IFLFRM(IGR2),
     &                    IFLIDL(IGR2), IFLFAC(IGR2), IFLNHD(IGR2), 0)
             IF (STPNOW()) RETURN                                         34.01
!ADC             ENDIF
          END IF
          CALL SWBROADC(ARR2,MXG(IGR2)*MYG(IGR2),SWREAL)                  40.31 40.30
        ENDIF                                                             40.03
      ENDIF
!     Interpolation over the computational grid
!     structured grid
      DO 230 IX = 1, MXC
        DO 240 IY = 1, MYC
          INDX = KGRPNT(IX,IY)
          IF (INDX.GT.1) THEN                                             40.00
            XP = XCGRID(IX,IY)
            YP = YCGRID(IX,IY)
            UU = SVALQI (XP, YP, IGR1, ARR, 0, IX, IY)                    40.31 30.90
            IF (IGR2.EQ.0) THEN
              COMPDA(INDX,JX3) = UU
            ELSE
              VV = SVALQI (XP, YP, IGR2, ARR2, 0, IX, IY)                 40.31 30.90
              COMPDA(INDX,JX3) =  UU*COSFC + VV*SINFC
              COMPDA(INDX,JY3) = -UU*SINFC + VV*COSFC
            ENDIF
          ENDIF
 240    CONTINUE
 230  CONTINUE
!     unstructured grid
      DO INDX = 1, nverts                                                 40.80
         XP = xcugrd(INDX)
         YP = ycugrd(INDX)
         IF ( IGTYPE(IGR1).EQ.3 ) THEN
            UU = ARR(INDX)
         ELSE
            UU = SVALQI (XP, YP, IGR1, ARR, 0, 0, 0)
         ENDIF
         IF (IGR2.EQ.0) THEN
            COMPDA(INDX,JX3) = UU
         ELSE
            IF ( IGTYPE(IGR2).EQ.3 ) THEN
               VV = ARR2(INDX)
            ELSE
               VV = SVALQI (XP, YP, IGR2, ARR2, 0, 0, 0)
            ENDIF
            COMPDA(INDX,JX3) =  UU*COSFC + VV*SINFC
            COMPDA(INDX,JY3) = -UU*SINFC + VV*COSFC
         ENDIF
      ENDDO                                                               40.80
      GOTO 200
!
!         Interpolation in time
!
 400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
      W1 = 1.-W3
      IF (ITEST.GE.60) WRITE(PRTEST,402) IGR1,
     &        TIMCO,IFLTIM(IGR1),W1,W3,JX1,JY1,JX2,JY2,JX3,JY3
 402  FORMAT (' input field', I2, ' interp at ', 2F9.0, 2F8.3, 6I3)
      DO 500 INDX = 1, MCGRD
        UU = W1 * COMPDA(INDX,JX2) + W3 * COMPDA(INDX,JX3)
        IF (IGR2.LE.0) THEN
          COMPDA(INDX,JX2) = UU
        ELSE
          VV = W1 * COMPDA(INDX,JY2) + W3 * COMPDA(INDX,JY3)
          VTOT = SQRT (UU*UU + VV*VV)
!
!         procedure to prevent loss of magnitude due to interpolation
!
          IF (VTOT.GT.0.) THEN
            SIZE1 = SQRT(COMPDA(INDX,JX2)**2 + COMPDA(INDX,JY2)**2)
            SIZE3 = SQRT(COMPDA(INDX,JX3)**2 + COMPDA(INDX,JY3)**2)
            SIZE2 = W1*SIZE1 + W3*SIZE3
!           SIZE2 is to be length of vector
            COMPDA(INDX,JX2) = SIZE2*UU/VTOT
            COMPDA(INDX,JY2) = SIZE2*VV/VTOT
          ELSE
            COMPDA(INDX,JX2) = UU
            COMPDA(INDX,JY2) = VV
          ENDIF
        ENDIF
 500  CONTINUE
      RETURN
!
!     End of subroutine FLFILE
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWINCO (AC2    ,COMPDA ,
     &                   XCGRID ,YCGRID ,                                 30.72
     &                   KGRPNT ,SPCDIR ,
     &                   SPCSIG ,XYTST   )                                30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72

      REAL    FDLSS,  TPDLSS, HSDLSS                                      30.80
!
      REAL     COMPDA(MCGRD,MCMVAR)
!
      REAL     AC2(MDC,MSC,MCGRD)
!
      INTEGER  KGRPNT(MXC,MYC), XYTST(*)                                  30.70
!
      LOGICAL  INTERN
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SWINCO')
!
!     *** Fetch Computation 'the mean delta' ***
      IF (KSPHER.EQ.0) THEN
        TLEN  = (XCLEN + YCLEN)/2.
      ELSE
        COSYG = COS (DEGRAD * (YOFFS + 0.5*(YCGMIN+YCGMAX)))              33.09
        TLEN  = LENDEG * (COSYG*XCLEN + YCLEN)/2.
      ENDIF
      TDXY  = FLOAT(MXCGL + MYCGL)/2.                                     40.31
      IF ( nvertsg.NE.0 ) TDXY = FLOAT(nvertsg)                           40.95 40.80
!
      FETCH = TLEN/TDXY
      SY0   = 3.3                                                         30.70
      IF (ITEST.GE.60 .OR. NPTST.GT.0) WRITE (PRTEST, 62) FETCH           30.70
  62  FORMAT (' test SWINCO, fetch:', E12.4)                              30.70
!
!     --- structured grid
!
      IF (ONED) THEN                                                      40.13
        IY1 = 1                                                           40.13
        IY2 = 1                                                           40.13
      ELSE                                                                40.13
        IY1 = 2                                                           40.13
        IY2 = MYC-1                                                       40.13
      ENDIF                                                               40.13
      DO IX = 2, MXC-1                                                    40.00
        DO IY = IY1, IY2                                                  40.13
!         check if the point is a true internal point                     40.00
          INTERN = .TRUE.
          INX = KGRPNT(IX-1,IY)                                           40.00
          IF (INX.LE.1) INTERN = .FALSE.                                  40.00
          INX = KGRPNT(IX+1,IY)                                           40.00
          IF (INX.LE.1) INTERN = .FALSE.                                  40.00
          IF (.NOT.ONED) THEN                                             40.13
            INX = KGRPNT(IX,IY-1)                                         40.00
            IF (INX.LE.1) INTERN = .FALSE.                                40.00
            INX = KGRPNT(IX,IY+1)                                         40.00
            IF (INX.LE.1) INTERN = .FALSE.                                40.00
          ENDIF                                                           40.13
          INX = KGRPNT(IX,IY)                                             30.70
          IF (INX.LE.1) INTERN = .FALSE.                                  40.00
          IF (INTERN) THEN                                                40.00
            TESTFL = .FALSE.                                              30.70
            DO IPTST = 1, NPTST                                           30.70
              IF (IX.EQ.XYTST(2*IPTST-1) .AND.                            30.70
     &            IY.EQ.XYTST(2*IPTST)) TESTFL = .TRUE.                   30.70
            ENDDO                                                         30.70
!
            IF (VARWI) THEN
              WX  = COMPDA(INX,JWX2)
              WY  = COMPDA(INX,JWY2)
!
!             *** Local wind speed and direction ***
              WSLOC = SQRT(WX*WX + WY*WY)
              IF (WX .NE. 0. .OR. WY .NE. 0.) THEN
                WDLOC = ATAN2(WY,WX)
              ELSE
                WDLOC = 0.
              ENDIF
            ELSE
!             uniform wind field
              WSLOC = U10                                                 30.60
              WDLOC = WDIP                                                30.60
            ENDIF
!
            IF (WSLOC .GT. 1.E-10) THEN                                   30.60
!
! Dimensionless Hs and Tp calculated according to K.K. Kahma & C.J. Calkoen,
! (JPO, 1992) and Pierson-Moskowitz for limit values.
!
!             calculate dimensionless fetch:
              FDLSS = GRAV * FETCH / (WSLOC*WSLOC)                        30.82
!
!             calculate dimensionless significant wave height:
              HSDLSS = MIN (0.21, 0.00288*FDLSS**0.45)                    40.00
              SPPARM(1) = HSDLSS * WSLOC**2 / GRAV                        40.00
!             calculate dimensionless peak period:
              TPDLSS = MIN (1./0.13, 0.46*FDLSS**0.27)                    40.00
              SPPARM(2) = WSLOC * TPDLSS / GRAV                           40.00
              IF (SPPARM(1).LT.0.05) SPPARM(2) = 2.                       40.51
              SPPARM(3) = 180. * WDLOC / PI
              SPPARM(4) = 2.
              IF (TESTFL) WRITE (PRTEST, 65) XCGRID(IX,IY),               30.70
     &        YCGRID(IX,IY), FDLSS, (SPPARM(JJ), JJ = 1, 3)               30.70
  65          FORMAT (' test point ',6(1X,E12.4))                         30.70
            ELSE
              SPPARM(1) = 0.02
              SPPARM(2) = 2.                                              40.51
              SPPARM(3) = 0.
              SPPARM(4) = 0.
            ENDIF
            CALL SSHAPE (AC2(1,1,INX), SPCSIG, SPCDIR, 2, 2)              40.00
!
          ENDIF
        ENDDO
      ENDDO
!
!     --- unstructured grid
!
      DO INX = 1, nverts                                                  40.80
!        internal vertices and ghost vertices only
         IF ( vmark(INX) == 0 .or. vmark(INX) == 999 ) THEN               40.95
            TESTFL = .FALSE.
            DO IPTST = 1, NPTST
              IF (INX.EQ.XYTST(IPTST)) TESTFL = .TRUE.
            ENDDO
!
            IF (VARWI) THEN
              WX = COMPDA(INX,JWX2)
              WY = COMPDA(INX,JWY2)
!
!             *** Local wind speed and direction ***
              WSLOC = SQRT(WX*WX + WY*WY)
              IF (WX .NE. 0. .OR. WY .NE. 0.) THEN
                WDLOC = ATAN2(WY,WX)
              ELSE
                WDLOC = 0.
              ENDIF
            ELSE
!             uniform wind field
              WSLOC = U10
              WDLOC = WDIP
            ENDIF
!
            IF (WSLOC .GT. 1.E-10) THEN
!
! Dimensionless Hs and Tp calculated according to K.K. Kahma & C.J. Calkoen,
! (JPO, 1992) and Pierson-Moskowitz for limit values.
!
!             calculate dimensionless fetch:
              FDLSS = GRAV * FETCH / (WSLOC*WSLOC)
!
!             calculate dimensionless significant wave height:
              HSDLSS = MIN (0.21, 0.00288*FDLSS**0.45)
              SPPARM(1) = HSDLSS * WSLOC**2 / GRAV
!             calculate dimensionless peak period:
              TPDLSS = MIN (1./0.13, 0.46*FDLSS**0.27)
              SPPARM(2) = WSLOC * TPDLSS / GRAV
              IF (SPPARM(1).LT.0.05) SPPARM(2) = 2.
              SPPARM(3) = 180. * WDLOC / PI
              SPPARM(4) = 2.
              IF (TESTFL) WRITE (PRTEST, 65) xcugrd(INX),
     &        ycugrd(INX), FDLSS, (SPPARM(JJ), JJ = 1, 3)
            ELSE
              SPPARM(1) = 0.02
              SPPARM(2) = 2.
              SPPARM(3) = 0.
              SPPARM(4) = 0.
            ENDIF
            CALL SSHAPE (AC2(1,1,INX), SPCSIG, SPCDIR, 2, 2)
!
         ENDIF
      ENDDO
!
      RETURN
! * end of subroutine SWINCO *
      END
!****************************************************************
!
      SUBROUTINE SWCLME
!
!****************************************************************
!
      USE M_WCAP
      USE M_GENARR
      USE M_PARALL
      USE M_DIFFR
      USE SwanGriddata
      USE SwanCompdata

      IMPLICIT NONE

      IF (ALLOCATED(SIGPOW))   DEALLOCATE(SIGPOW)
      IF (ALLOCATED(KGRPNT))   DEALLOCATE(KGRPNT)
      IF (ALLOCATED(KGRBND))   DEALLOCATE(KGRBND)
      IF (ALLOCATED(XYTST ))   DEALLOCATE(XYTST )
      IF (ALLOCATED(AC2   ))   DEALLOCATE(AC2   )
      IF (ALLOCATED(XCGRID))   DEALLOCATE(XCGRID)
      IF (ALLOCATED(YCGRID))   DEALLOCATE(YCGRID)
      IF (ALLOCATED(SPCSIG))   DEALLOCATE(SPCSIG)
      IF (ALLOCATED(SPCDIR))   DEALLOCATE(SPCDIR)
      IF (ALLOCATED(DEPTH ))   DEALLOCATE(DEPTH )
      IF (ALLOCATED(FRIC  ))   DEALLOCATE(FRIC  )
      IF (ALLOCATED(UXB   ))   DEALLOCATE(UXB   )
      IF (ALLOCATED(UYB   ))   DEALLOCATE(UYB   )
      IF (ALLOCATED(WXI   ))   DEALLOCATE(WXI   )
      IF (ALLOCATED(WYI   ))   DEALLOCATE(WYI   )
      IF (ALLOCATED(WLEVL ))   DEALLOCATE(WLEVL )
      IF (ALLOCATED(ASTDF ))   DEALLOCATE(ASTDF )
      IF (ALLOCATED(IBLKAD))   DEALLOCATE(IBLKAD)
      IF (ALLOCATED(XGRDGL))   DEALLOCATE(XGRDGL)
      IF (ALLOCATED(YGRDGL))   DEALLOCATE(YGRDGL)
      IF (ALLOCATED(KGRPGL))   DEALLOCATE(KGRPGL)
      IF (ALLOCATED(KGRBGL))   DEALLOCATE(KGRBGL)
      IF (ALLOCATED(DIFPARAM)) DEALLOCATE(DIFPARAM)
      IF (ALLOCATED(DIFPARDX)) DEALLOCATE(DIFPARDX)
      IF (ALLOCATED(DIFPARDY)) DEALLOCATE(DIFPARDY)
      IF (ALLOCATED(NPLAF ))   DEALLOCATE(NPLAF )
      IF (ALLOCATED(LAYH  ))   DEALLOCATE(LAYH  )
      IF (ALLOCATED(VEGDIL))   DEALLOCATE(VEGDIL)
      IF (ALLOCATED(VEGNSL))   DEALLOCATE(VEGNSL)
      IF (ALLOCATED(VEGDRL))   DEALLOCATE(VEGDRL)
!
      IF (ALLOCATED(xcugrd))   DEALLOCATE(xcugrd)
      IF (ALLOCATED(ycugrd))   DEALLOCATE(ycugrd)
      IF (ALLOCATED( vmark))   DEALLOCATE( vmark)
      IF (ALLOCATED( vlist))   DEALLOCATE( vlist)
      IF (ALLOCATED( blist))   DEALLOCATE( blist)
!
      RETURN
      END
