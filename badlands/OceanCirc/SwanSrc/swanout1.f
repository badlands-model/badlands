!
!     SWAN/OUTPUT       file 1 of 2
!
!  Contents of this file:
!     SWOUTP
!     SWORDC
!     SWODDC
!     SWOEXC
!     SWOEXD
!     SWIPOL
!     SWOEXA
!     SWOINA
!     SWOEXF
!
!     main output routine and computation of output quantities
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOUTP ( AC2             ,                                40.31 30.90
     &                   SPCSIG          ,SPCDIR  ,                       30.72
     &                   COMPDA          ,XYTST   ,
     &                   KGRPNT          ,XCGRID  ,                       30.72
     &                   YCGRID          ,OURQT , FIELD  )                       40.51 40.30
!                                                                      *
!***********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.80
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
      USE swan_coupler_export
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.81: Annette  Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     34.01: Jeroen Adema
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.86: Nico Booij
!     40.90: Nico Booij
!
!  1. Updates
!
!     10.10, Aug. 94: computation of force is added (subr. SWOEXF)
!                     arrays NE and NED added
!     30.72, Oct. 97: changed floating point comparison to avoid equality
!                     comparisons
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, June 98: argument KGRBND added, call SWPLOT and SWOEXC
!                     modified
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Jan. 99: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     40.00, Jan. 99: argument RTYPE added in call SWODDC
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Oct. 00: Made TYPE of several equivalenced arrays correct
!     40.02, Oct. 00: Modified argument list of SWPLOT to avoid int/real conflict
!     40.13, Oct. 01: Forces always computed by post-processing procedure
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing POOL construction and HPGL-functionality
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: optimization output process in parallel mode
!     40.80, Feb. 08: computation of wave-induced force on unstructured grid added
!     40.86, Feb. 08: arguments added to calls of subroutines
!                     to prevent interpolation over obstacles
!     40.90, June 08: arguments added to call of subroutine SWSPEC
!                     to prevent interpolation over obstacles
!
!  2. Purpose
!
!     Processing of the output requests
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! i   OURQT : array indicating at what time requested output              40.51
!             is processed                                                40.51
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
!
      REAL    OURQT(MAX_OUTP_REQ)                                         40.51 40.30
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
      REAL    FIELD(2,MXCGL*MYCGL)
!
!     AC2     real arr input  action density in all computational points
!     SPCDIR  real arr input  spectral directions, cosines and sines
!     UX2     real arr input  current velocity x-comp.
!     UY2     real arr input  current velocity y-comp.
!     UBOT    real arr input  orbital velocities at bottom
!     WX2     real arr input  wind velocity x-comp.
!     WY2     real arr input  wind velocity y-comp.
!
!  8. Subroutines used
!
!     SWBLOK
!     SWTABP
!     SWSPEC
!     FOR
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     MAIN program
!
! 10. Error messages
!
!     If an output request is of an unknown type, an error message
!     is printed.
!
! 11. Remarks
!
!     Data:
!
!     output requests are encoded in array OUTREQ; these are set by
!     commands TABLE, BLOCK, SPEC, etc. (see subr SWREOQ in file SWANPRE2)
!     each output request refers to one set of output locations,
!     and to one or more output quantities.
!     1st value in OUTREQ: time of next output, 2nd value: interval between
!     outputs, 3d value: type of output request RTYPE (encoded as integer),
!     4&5: PSNAME (name of point set),
!     6: file unit number, 7..10: output filename,
!     other: dependent on type of output.
!
!     data on output locations are in array OUTDA; these are set by
!     commands FRAME, POINTS, CURVE etc. (see subr SWREPS in file SWANPRE2)
!     each set is characterized by its name (SNAME in the code)
!     STYPE is the type of set (i.e. 'F' for Frame etc.)
!
!     properties of output quantities are in arrays OVSNAM, OVLNAM,
!     OVUNIT, OVSVTY etc.; these are set in subr SWINIT (file SWANMAIN)
!     each output quantity is assigned a fixed number; i.e. 1=Xp, 2=Yp,
!     7=Dissip, 10=Hs, 11=Tm01 etc.
!     subr SVARTP determines the above number from the name of the
!     quantity as it appears in the user command; this is compared with
!     OVKEYW.
!
!     Procedure:
!
!     After the coordinates of all output locations have been determined,
!     values of all output quantities are calculated, and written into
!     2d array VOQ (one or two columns for each output quantity, one line
!     for each location). array SWOUTP shows with quantity is written in
!     each column.
!     After array VOQ is filled, the actual output starts; which subroutine
!     is called depends on RTYPE (see structure scheme below).
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     For all output requests do
!         Call SWORDC (analyzes output req.; gives name of output
!                      point set, type of output RTYPE
!                      and number of quantities per point)
!         Call SWODDC (analyzes output data; gives number of output
!                      points)
!         ------------------------------------------------------------
!         Call SWOEXC (compute coordinates of output points)
!         Call SWOEXD (compute depth, current vel. dissipation etc.)
!         Call SWOEXA (compute action density and related quant.)
!         Call SWOEXF (compute wave-driven force)
!         ------------------------------------------------------------
!         If RTYPE = 'BLKP', 'BLKD' or 'BLKL' then
!                     call SWBLOK for block output
!         If RTYPE = 'TABP' or 'TABD' then call SWTABP for output in
!                     table
!         If RTYPE = 'SPEC' then call SWSPEC for spectral output
!     ----------------------------------------------------------------
!     If program is run in stationary mode                                40.00
!     Then Close all opened files                                         40.00
!     ----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER   VOQR(NMOVAR)   ,BKC           ,                           40.31
     &          XYTST(*)       ,KGRPNT(MXC,MYC),IERR                      30.81 30.21
!
      REAL      AC2(MDC,MSC,MCGRD) ,
     &          COMPDA(MCGRD,MCMVAR)
      LOGICAL, ALLOCATABLE :: CROSS(:,:) ! true if obstacle is between    40.86
                                         ! output point and computational 40.86
                                         ! grid point                     40.86
!
      INTEGER, ALLOCATABLE :: IONOD(:)                                    40.51
      REAL, ALLOCATABLE :: ACLOC(:), AUX1(:), VOQ(:)                      40.31
      REAL, ALLOCATABLE :: FORCE(:,:)                                     40.80
      REAL              :: KNUM(MSC), CG(MSC), NE(MSC), NED(MSC)          40.31
      TYPE(OPSDAT), POINTER :: CUOPS                                      40.31
      TYPE(ORQDAT), POINTER :: CORQ                                       40.31
!
      LOGICAL   OQPROC(NMOVAR)  , LOGACT                                  30.00
      CHARACTER RTYPE *4, STYPE *1, PNAME *8, PTYPE *1
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWOUTP')
!
!     processing of output requests
!
      IF (NPTST.GT.0) CALL AC2TST (XYTST, AC2 ,KGRPNT)                    30.21
!
      IF (NREOQ.EQ.0) THEN
        CALL MSGERR (1, 'no output requested')                            40.31
        RETURN
      ENDIF
      IF (ITEST.GE.10) WRITE (PRINTF, 12) NREOQ
  12  FORMAT (1X, I3, ' output requests')
!
!     Repeat for all output requests:
!
      CORQ => FORQ                                                        40.31
      DO 70 IRQ = 1, NREOQ
!
!       ***** processing of output instructions *****
!
        BKC   = 0
        NVOQP = 0
!
!       call SWORDC to analyse output request encoded in array OUTREQ
!       NVOQP (number of output quantities)
!
        RTYPE = CORQ%RQTYPE                                               40.31
        SNAME = CORQ%PSNAME                                               40.31
        CALL SWORDC (CORQ%OQI, CORQ%OQR, CORQ%IVTYP, RTYPE,               40.31 30.90
     &               SNAME, NVOQP, OQPROC, BKC, VOQR, OURQT(IRQ),         40.51
     &               LOGACT)                                              30.00
        IF (.NOT.LOGACT) THEN                                             40.31 30.00
           CORQ => CORQ%NEXTORQ                                           40.31
           GOTO 70                                                        40.31 30.00
        END IF
!
!        IF (SCREEN.NE.PRINTF.AND.IAMMASTER) WRITE (SCREEN, 15) IRQ        40.30
!  15    FORMAT ('+SWAN is processing output request ', I4)                40.41 30.00
        IF (ITEST.GE.10) WRITE(PRINTF,16) IRQ
  16    FORMAT (' SWAN is processing output request ', I4)                40.41
!
        CUOPS => FOPS                                                     40.31
        DO                                                                40.31
          IF (CUOPS%PSNAME.EQ.SNAME) EXIT                                 40.31
          IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) THEN                        40.31
             CALL MSGERR (3, 'Output requested for non-existing points')  10.21
             WRITE (PRINTF, 18) SNAME                                     40.31 30.81 30.00
  18         FORMAT (' Point set: ', A)                                   40.31 30.00
             GOTO 68                                                      40.00
          END IF                                                          40.31
          CUOPS => CUOPS%NEXTOPS                                          40.31
        END DO                                                            40.31
!
        IF (ITEST.GE.80 .OR. IOUTES .GE. 10)
     &  WRITE (PRTEST, 22) IRQ, NVOQP, RTYPE, SNAME
  22    FORMAT (' Test SWOUTP ', 2I6, 2X, A4, 2X, A16)
!
!       call SWODDC to analyse output data; results: STYPE (type of output
!       point set), MIP (number of output locations) etc.
!
        STYPE = CUOPS%PSTYPE                                              40.31
        MIP   = CUOPS%MIP                                                 40.31
        CALL SWODDC (CUOPS%OPI, CUOPS%OPR, SNAME, STYPE, MIP, MXK,        40.31
     &               MYK, XNLEN, YNLEN, MXN, MYN, XPCN, YPCN, ALPCN,
     &               XCGRID,YCGRID,RTYPE)                                 40.00
!
!       assign memory to array VOQ (contains output quantities for all
!                                   output points)
        ALLOCATE(VOQ(MIP*NVOQP))                                          40.31
!
!       assign memory to array CROSS (indicates crossing of obstacles     40.86
!                                     in between output and grid points)  40.86
        ALLOCATE(CROSS(4,MIP))                                            40.86
!
!       assign memory to array IONOD (indicates in which subdomain        40.51
!                                     output points are located)          40.51
        ALLOCATE(IONOD(MIP))                                              40.51
        IONOD = -999                                                      40.51
!
!       call SWOEXC to calculate quantities dependent only on coordinates
!
        CALL SWOEXC (STYPE               ,                                40.31
     &               CUOPS%OPI           ,CUOPS%OPR           ,           40.31
     &               CUOPS%XP            ,CUOPS%YP            ,           40.31
     &               MIP                 ,VOQ(1)              ,           40.31 30.90
     &               VOQ(1+MIP)          ,VOQ(1+2*MIP)        ,           40.31 30.90
     &               VOQ(1+3*MIP)        ,KGRPNT              ,           40.31 30.90
     &               XCGRID              ,YCGRID              ,           30.21
     &               CROSS                                    )           40.86 40.00
!
!       Compute wave-induced force on unstructured grid                   40.80
!
        IF (OQPROC(20) .AND. OPTG.EQ.5) THEN                              40.80
           ALLOCATE(FORCE(nverts,2))                                      40.80
           CALL SwanComputeForce ( FORCE(1,1), FORCE(1,2), AC2,           40.80
     &                             COMPDA(1,JDP2), SPCSIG, SPCDIR )       40.80
        ELSE                                                              40.80
           ALLOCATE(FORCE(0,0))                                           40.80
        ENDIF                                                             40.80
!
!       call SWOEXD to interpolate quantities which are computed during the
!       SWAN computation, such as Qb, Dissipation, Ursell etc.
!
        CALL SWOEXD (OQPROC, MIP, VOQ(1+2*MIP),                           40.31 30.90
     &               VOQ(1+3*MIP), VOQR, VOQ(1),                          40.31 30.90
     &               COMPDA, KGRPNT, FORCE, CROSS, IONOD)                 40.86 40.80 40.31
!
        DEALLOCATE(FORCE)                                                 40.80
!
        IF (BKC .GT. 0) THEN
!
!         assign memory to array ACLOC (contains spectrum for one output point)
!
          ALLOCATE(ACLOC(MDC*MSC))                                        40.31
!
!         call SWOEXA to compute quantities for which spectrum is needed
!         (except wave-induced force)
!
          CALL SWOEXA (OQPROC              ,BKC                 ,
     &                 MIP                 ,VOQ(1+2*MIP)        ,         40.31 30.90
     &                 VOQ(1+3*MIP)        ,VOQR                ,         40.31 30.90
     &                 VOQ(1)              ,AC2                 ,         40.31 30.90
     &                 ACLOC               ,SPCSIG              ,         40.31 30.90
     &                 KNUM                ,CG                  ,         40.31 30.90
     &                 SPCDIR              ,NE                  ,         40.31 30.90
     &                 NED                 ,KGRPNT              ,         40.31 30.90
     &                 COMPDA(1,JDP2)      ,CROSS               )         40.86
!
!         call SWOEXF to compute wave-driven force on regular grid
!
          IF (OQPROC(20) .AND. OPTG.NE.5)                                 40.80 40.13
     &      CALL SWOEXF (MIP                 ,VOQ(1+2*MIP)         ,      40.31 30.90
     &                   VOQ(1+3*MIP)        ,VOQR                 ,      40.31 30.90
     &                   VOQ(1)              ,AC2                  ,      40.31 30.90
     &                   COMPDA(1,JDP2)      ,SPCSIG               ,      30.72
     &                   KNUM                ,CG                   ,      40.31 30.90
     &                   SPCDIR              ,NE                   ,      40.31 30.90
     &                   NED                 ,KGRPNT               ,      40.31 30.90
     &                   XCGRID              ,YCGRID               ,      30.72
     &                   IONOD                                            40.31
     &                                                             )
!
          DEALLOCATE(ACLOC)                                               40.31
        ENDIF
!
        IF (ITEST.GE.100 ) THEN
          WRITE (PRTEST, 23) (VOQR(II), II=1, NMOVAR)
  23      FORMAT (' arrays VOQR and VOQ:', 30I3)
          DO 25 IP=1, MIN(MIP,20)
            WRITE (PRTEST, 24) (VOQ(IP+(JJ-1)*MIP),
     &                          JJ=1, NVOQP)
  24        FORMAT (12(1X,E10.4))
  25      CONTINUE
        ENDIF
!
!       ***** block output *****
        IF (RTYPE(1:3) .EQ. 'BLK') THEN
          IF (PARLL) THEN                                                 40.31
             CALL SWBLKP ( CORQ%OQI, CORQ%IVTYP, MXK, MYK, VOQR,          40.31
     &                     VOQ(1), IONOD )                                40.51 40.31
          ELSE                                                            40.31
             CALL SWBLOK ( RTYPE, CORQ%OQI, CORQ%IVTYP, CORQ%FAC,         40.31
     &                     SNAME, MXK, MYK, IRQ, VOQR, VOQ(1) )           40.51 40.31
          END IF                                                          40.31
          IF (STPNOW()) RETURN                                            34.01
          GOTO 68                                                         40.00
        ENDIF
!
!       ***** table output *****
        IF (RTYPE(1:3) .EQ. 'TAB') THEN
          CALL SWTABP ( RTYPE, CORQ%OQI, CORQ%IVTYP, SNAME,               40.31
     &                  MIP, VOQR, VOQ(1), IONOD )                        40.51 40.31
          IF (STPNOW()) RETURN                                            34.01
          GOTO 68                                                         40.00
        ENDIF
!
!       ***** spectral output *****
        IF (RTYPE(1:2) .EQ. 'SP') THEN                                    20.28
          IF (RTYPE(4:4).EQ.'C') THEN
             ALLOCATE(AUX1(MSC*MDC))                                      40.31
          ELSE
             ALLOCATE(AUX1(3*MSC))                                        40.31
          ENDIF
          CALL SWSPEC ( RTYPE, CORQ%OQI, MIP, VOQR, VOQ(1), AC2, AUX1,    40.31
     &                  SPCSIG, SPCDIR, COMPDA(1,JDP2), KGRPNT, CROSS,    40.90 40.31
     &                  IONOD )                                           40.31
          IF (STPNOW()) RETURN                                            34.01
          DEALLOCATE(AUX1)                                                40.31
          GOTO 68                                                         40.00
        ENDIF
!
  60    WRITE (PRINTF, 62) IRQ, NVOQP, RTYPE, SNAME, MIP
  62    FORMAT (' Error in output request ', 2I6, 2X, A4, 2X, A16, I6)
!
  68    CONTINUE                                                          40.31
 
        call export_swan(MIP,MXK,MYK,NVOQP,VOQ,VOQR,IONOD)

	    IF (STPNOW()) RETURN
	    
        DEALLOCATE(VOQ,CROSS,IONOD)                                       40.86 40.51 40.31
        CORQ => CORQ%NEXTORQ                                              40.31
  70  CONTINUE
!
!     Termination of output
!
  200 RETURN
      END SUBROUTINE SWOUTP

!***********************************************************************
!                                                                      *
      SUBROUTINE SWORDC (OUTI, OUTR, IVTYP, RTYPE, PSNAME, NVOQP,         40.31
     &                   OQPROC, BKC,                                     40.31
     &                   VOQR, OURQT, LOGACT)                             40.51 30.00
!                                                                      *
!***********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.61: Roberto Padilla
!     30.74: IJsbrand Haagsma (Include version)
!     30.81: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Update
!
!     10.33, Jan. 95: computation of wind added in polar plot
!     30.61, Summ 97: give values for array VOQR when OQPROC = .TRUE.
!                     in case of 'nest'
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.81, Jan. 99: Replaced variable FROM by FROM_ (because FROM
!                     is a reserved word)
!     40.30, May  03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing HPGL-functionality
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Decodes output requests
!
!  3. Method
!
!     ---
!
!  4. Argument list
!
!     OUTR    Real   input    Code for one output request
!     RTYPE   Char   outp     type of output
!     PSNAME  Char   outp     name of output point set referred to
!     NVOQP   Int    outp     number of data per output point
!     OQPROC  logic  outp     whether or not an output quantity  must
!                             be processed
!     VOQR    Int ar outp     place of each output quantity
!                             (subscript: IVTYP)
!     OURQT   Int ar input    array indicating at what time requested     40.51
!                             output is processed                         40.51
!
!  8. Subroutines used
!
!     SPSET
!     SUVIPL
!     SBLKPT
!     SCUNIT
!     SFLFUN (all SWAN/OUTP)
!     TABHED
!     MSGERR
!     FOR
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     If the point set is not of the type frame, an error message
!     is printed and control returns to subroutine OUTPUT
!
! 11. Remarks
!
!     output interval negative means that output is made only at end
!     of computation                                                      40.00
!
! 12. Structure
!
!     -----------------------------------------------------------------   40.00
!     If dynamic mode
!     Then determine TNEXT (time of next requested output)
!          determine DIF (interval between end time and present time)
!          If DIF is less than half time step and output interval is
!               negative
!          Then enable output (by making LOGACT = true)
!          Else If time of computation >= TNEXT and output interval is
!               positive
!          Then enable output
!          Else disable output (by making LOGACT = false)
!               Return
!          ------------------------------------------------------------   40.00
!     Else enable output
!     -----------------------------------------------------------------   40.00
!     Set all OQPROC = false (if OQPROC is true corresponding quantity
!                             must be computed)
!     Make OQPROC true for quantities Xp, Yp, Xc and Yc
!     Set all values of VOQR = 0 (VOQR indicates where value of a
!                                 quantity is stored in array VOQ)
!     Make VOQR nonzero for quantities Xp, Yp, Xc and Yc
!     Assign value to NVAR depending on type of output request
!     -----------------------------------------------------------------   40.00
!
! 13. Source text
!
      INTEGER    VOQR(*), OUTI(*), BKC                                    30.00
      INTEGER    IVTYP(*)                                                 40.31
      REAL       OUTR(*)                                                  30.00
      REAL       OURQT                                                    40.51
      LOGICAL    OQPROC(NMOVAR), LOGACT                                   30.00
      CHARACTER  PSNAME *(*), RTYPE *(*)                                  40.31 30.81
      SAVE IENT
      DATA IENT /0/                                                       40.31 30.81
      CALL STRACE (IENT, 'SWORDC')
!
!     check time of output action:                                        30.00
      IF (NSTATM.EQ.1) THEN                                                    40.00
!       check time of output action:                                      40.00
!       DIF  in case that timco is not a fraction of the
!       computational period and the user do not ask for a periodic plots
        DIF = TFINC - TIMCO
        TNEXT = OUTR(1)
        IF ( PARLL.AND.OURQT.EQ.-9999.) OURQT = OUTR(1)                   40.51 40.30
        IF (ITEST.GE.60) WRITE (PRTEST, *) ' output times ', TNEXT,
     &        OUTR(2), DT, TFINC, TIMCO
        IF (ABS(DIF).LT.0.5*DT .AND. OUTR(2).LT.0.) THEN                  40.00
          OUTR(1) = TIMCO
          LOGACT = .TRUE.
        ELSE IF (OUTR(2).GT.0. .AND. TIMCO.GE.TNEXT) THEN                 40.00
          OUTR(1) = TNEXT + OUTR(2)                                       40.00
          LOGACT = .TRUE.
        ELSE
          LOGACT = .FALSE.
          RETURN
        ENDIF
      ELSE
        LOGACT = .TRUE.                                                   30.00
      ENDIF                                                               30.00
!
!     action is taken, proceed with analysing output request
!
!     Ivtype 1 and 2 are Xp and Yp
!
      OQPROC(1) = .TRUE.
      VOQR(1)   = 1
      OQPROC(2) = .TRUE.
      VOQR(2)   = 2
!
!     clear VOQR and OQPROC from old information
!
      DO 10 IVT = 3, NMOVAR
        VOQR(IVT)   = 0
        OQPROC(IVT) = .FALSE.
  10  CONTINUE
!
!     Ivtype 24 and 25 are Xc and Yc
!
      OQPROC(24) = .TRUE.
      VOQR(24)   = 3
      OQPROC(25) = .TRUE.
      VOQR(25)   = 4
      NVOQP = 4
!
      NVAR = OUTI(3)                                                      40.31
!
      DO 20 IVAR = 1, NVAR
         IVTYPE = IVTYP(IVAR)                                             40.31
!
         IF (IVTYPE.LT.1 .OR. IVTYPE.GT.NMOVAR) THEN
           CALL MSGERR (2, 'wrong value for IVTYPE')
           WRITE (PRINTF, 17) RTYPE, PSNAME, IVTYPE, NVAR
  17       FORMAT (' type, points, var: ', A4, 2X, A8, 2X, 2I8)
           GOTO 20
         ENDIF
!
         IF (OVSVTY(IVTYPE).LE.2 .AND. .NOT.OQPROC(IVTYPE)) THEN
!           output quantity is a scalar
            NVOQP = NVOQP + 1
            VOQR(IVTYPE) = NVOQP
         ELSE IF (OVSVTY(IVTYPE).EQ.3 .AND. .NOT.OQPROC(IVTYPE)) THEN
!           output quantity is a vector
            NVOQP = NVOQP + 2
            VOQR(IVTYPE) = NVOQP-1
         ENDIF
         OQPROC(IVTYPE) = .TRUE.
         IF (ITEST.GE.80 .OR. IOUTES .GE. 20) WRITE (PRTEST, 22)
     &   IVAR, IVTYPE, VOQR(IVTYPE)                                       10.09
  22     FORMAT (' SWORDC, output quantity:', 3I6)
!
!        for spectral width add Tm02 as output quantity                   20.61
!
         IF (IVTYPE.EQ.33) THEN
           IF (.NOT.OQPROC(32)) THEN
             NVOQP = NVOQP + 1
             VOQR(32) = NVOQP
             OQPROC(32) = .TRUE.
           ENDIF
         ENDIF
!
!        for BFI add steepness and Qp as output quantities                40.64
!
         IF (IVTYPE.EQ.59) THEN
           IF (.NOT.OQPROC(18)) THEN
             NVOQP = NVOQP + 1
             VOQR(18) = NVOQP
             OQPROC(18) = .TRUE.
           ENDIF
           IF (.NOT.OQPROC(58)) THEN
             NVOQP = NVOQP + 1
             VOQR(58) = NVOQP
             OQPROC(58) = .TRUE.
           ENDIF
         ENDIF
!
!        for some quantities compute action densities
!
         IF (IVTYPE.EQ.10 .OR. IVTYPE.EQ.11 .OR. IVTYPE.EQ.12 .OR.
     &       IVTYPE.EQ.13 .OR. IVTYPE.EQ.14 .OR. IVTYPE.EQ.16 .OR.
     &       IVTYPE.EQ.21 .OR. IVTYPE.EQ.22 .OR. IVTYPE.EQ.43 .OR.
     &       IVTYPE.EQ.44 .OR. IVTYPE.EQ.48 .OR. IVTYPE.EQ.53 .OR.
     &       IVTYPE.EQ.58                                          )      40.64 40.61 40.51 40.41 40.00
     &   BKC = MAX (1, BKC)
!
!        for some quantities also compute Depth, current, K and Cg
!
         IF (IVTYPE.EQ.15 .OR. IVTYPE.EQ.17 .OR. IVTYPE.EQ.18 .OR.
     &       IVTYPE.EQ.19 .OR. IVTYPE.EQ.20 .OR. IVTYPE.EQ.28 .OR.        10.10
     &       IVTYPE.EQ.32 .OR. IVTYPE.EQ.33 .OR. IVTYPE.EQ.42 .OR.
     &       IVTYPE.EQ.47 .OR. IVTYPE.EQ.59 .OR. IVTYPE.EQ.71      )      41.15 40.64 40.41 40.00
     &   BKC = 2
!
         IF (IVTYPE.EQ.11 .AND. ICUR.GT.0) BKC = 2                        20.36
         IF (BKC.GT.0) THEN
!           depth must be computed
            IF (.NOT.OQPROC(4)) THEN
               NVOQP = NVOQP + 1
               VOQR(4) = NVOQP
               OQPROC(4)=.TRUE.
            ENDIF
!           current velocity must be computed
            IF (.NOT.OQPROC(5) .AND. ICUR.GT.0) THEN
               NVOQP = NVOQP + 2
               VOQR(5) = NVOQP-1
               OQPROC(5)=.TRUE.
            ENDIF
         ENDIF
  20  CONTINUE
!
!     in case of print of spectrum Ux and Uy have to be computed          20.28
!
      IF (     RTYPE(1:2) .EQ. 'SP'                                       40.31
     &    .OR. RTYPE(1:2) .EQ. 'NE') THEN                                 30.61
        OQPROC(4)  = .TRUE.
        VOQR(4)    = NVOQP+1
        NVOQP      = NVOQP+1
        BKC        = 1
        IF (ICUR.GT.0) THEN
          BKC = 2
          OQPROC(5) = .TRUE.
          VOQR(5)   = NVOQP+1
          NVOQP     = NVOQP+2
        ENDIF
      ENDIF                                                               20.28
!
      RETURN
!*    end of subroutine SWORDC   **
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWODDC (OPI, OPR, PSNAME, PSTYPE, MIP, MXK, MYK,         40.31
     &                   XNLEN, YNLEN, MXN, MYN, XPCN, YPCN, ALPCN,
     &                   XCGRID,YCGRID,RTYPE)                             40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!            Oct. 95: New subroutine
!     30.72, Sept 97: Replaced DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, Jan. 99: argument RTYPE added, computation of ALCQ changed if
!                     output type (indicated by RTYPE) is PLOT
!     40.22, Sep. 01: small corrections
!     40.31, Nov. 03: removing HPGL-functionality
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Decodes output point set data
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!  PSNAME      Char   input   name of output point set referred to
!  PSTYPE      Char   input   type of output point set
!  MIP         Int    outp    number of output points
!  MXK         Int    outp    number of output points in X-direction (Frame)
!  MYK         Int    outp    number of output points in Y-direction (Frame)
!  XNLEN,YNLEN real   outp    (X,Y)lenght of the nested grid
!  MXN, MYN    int    outp    number of meshes in X, Y direction for
!                             the nested grid
!  XPCN, YPCN  real   outp    location of the origin of the nested grid
!  ALPCN       real   outp    angle of the nested grid with the positive
!                             x-axis, counterclockwise measured
!  RTYPE       char   input   indicates type of output; "PLOT" means that
!                             a spatial plot is made
      CHARACTER  RTYPE *(*)
      INTEGER OPI(2)                                                      40.31
      REAL    OPR(5)                                                      40.31
!
!  5. SUBROUTINES CALLING
!
!       SWOUTP (SWAN/OUTP)
!
!  6. SUBROUTINES USED
!
!  7. ERROR MESSAGES
!
!       If the point set is not of a known type an error message          40.00
!       is printed and control returns to subroutine SWOUTP
!       If the point set is not of the type frame or ngrid an error message
!       is printed and control returns to subroutine SWOUTP
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Depending on type of set of output points                         40.00
!       determine name, type and number of output points of
!       the output point set
!       ----------------------------------------------------------------
!
      CHARACTER  PSNAME *(*), PSTYPE *1
      LOGICAL    EQREAL
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWODDC')
!
      ALPQ  = 0.
      COSPQ = 1.
      SINPQ = 0.
!     ALCQ  = ALPC           removed 30.50: U and V now in user coordinates
      ALCQ  = 0.
      COSCQ = COS(ALCQ)
      SINCQ = SIN(ALCQ)
!
      IF (PSTYPE.EQ.'F') THEN
        MXK  = OPI(1)                                                     40.31
        MYK  = OPI(2)                                                     40.31
        MIP  = MXK * MYK
        XPQ  = OPR(1)                                                     40.31
        YPQ  = OPR(2)                                                     40.31
        ALPQ = OPR(5)                                                     40.31
        COSPQ = COS(ALPQ)
        SINPQ = SIN(ALPQ)
        XQP   = -XPQ*COSPQ - YPQ*SINPQ
        YQP   =  XPQ*SINPQ - YPQ*COSPQ
        IF (EQREAL(OUTPAR(4),1.)) THEN                                    40.00
!         directions will be w.r.t. frame coordinate system
          ALCQ = -ALPQ
        ELSE
!         directions will be w.r.t. user coordinate system (default)      40.00
          ALCQ = 0.                                                       40.00
        ENDIF                                                             40.00
        COSCQ = COS(ALCQ)
        SINCQ = SIN(ALCQ)
        DXK   = OPR(3) / FLOAT(MXK-1)                                     40.31
        IF ( MYK.GT.1 ) THEN                                              40.22
           DYK   = OPR(4) / FLOAT(MYK-1)                                  40.31 40.22
        ELSE                                                              40.22
           DYK   = 0.                                                     40.22
        END IF                                                            40.22
      ELSE IF (PSTYPE.EQ.'H') THEN                                        30.21
        MXK  = OPI(1)                                                     40.31
        MYK  = OPI(2)                                                     40.31
        MIP  = MXK * MYK
        ALPQ = OPR(5)                                                     40.31
        COSPQ = COS(ALPQ)
        SINPQ = SIN(ALPQ)
        XCMAX = OPR(1)                                                    40.31
        YCMAX = OPR(2)                                                    40.31
        XCMIN = OPR(3)                                                    40.31
        YCMIN = OPR(4)                                                    40.31
!       *** Find XQLEN and YQLEN taken the extreme points    ***
!       *** that belongs to the frame                        ***
        XPMIN =  1.E09
        YPMIN =  1.E09
        XPMAX = -1.E09
        YPMAX = -1.E09
!
        XPQ   = XPMIN
        YPQ   = YPMIN
        XQP   = -XPQ
        YQP   = -YPQ
        ALCQ  = 0.
        COSCQ = COS(ALCQ)
        SINCQ = SIN(ALCQ)
        DXK   = (XPMAX - XPMIN)/ FLOAT(MXK-1)
        IF ( MYK.GT.1 ) THEN                                              40.22
           DYK   = (YPMAX - YPMIN) / FLOAT(MYK-1)                         40.22
        ELSE                                                              40.22
           DYK   = 0.                                                     40.22
        END IF                                                            40.22
        XNLEN = 0.
        YNLEN = 0.
        MXN   = 0
        MYN   = 0
        XPCN  = 0.
        YPCN  = 0.
        ALPCN = 0.
        IF (IOUTES .GE. 20) THEN
          WRITE(PRINTF, 11) MXK ,MYK , XQP ,YQP , DXK ,DYK
 11       FORMAT ('SWODDC :',/,'     MXK ,MYK , XQP     ,YQP       ,DXK'
     &          ,'        ,DYK'
     &          ,/,2(1X,I5), 4(1X,E9.3))
          IF (PSTYPE .EQ. 'H') WRITE(PRINTF, 12)XCMIN,XCMAX,YCMIN,YCMAX,
     &      XPMAX ,XPMIN
 12         FORMAT(' XCMIN   ,XCMAX    ,YCMIN    ,YCMAX    ,',
     &      'XPMAX    ,XPMIN :',/,6(1X,E9.3))
        ENDIF
      ELSE IF (PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P') THEN
        MXK = 0
        MYK = 0
        XNLEN = 0.
        YNLEN = 0.
        MXN   = 0
        MYN   = 0
        XPCN  = 0.
        YPCN  = 0.
        ALPCN = 0.
      ELSE IF (PSTYPE.EQ.'N') THEN
!       nested grid                                                       40.00
        MXK = 0
        MYK = 0
        XNLEN = OPR(1)                                                    40.31
        YNLEN = OPR(2)                                                    40.31
        XPCN  = OPR(3)                                                    40.31
        YPCN  = OPR(4)                                                    40.31
        ALPCN = OPR(5)                                                    40.31
        MXN   = OPI(1)                                                    40.31
        MYN   = OPI(2)                                                    40.31
      ELSE IF (PSTYPE.EQ.'U') THEN                                        40.80
        MXK   = MIP                                                       40.80
        MYK   = 1                                                         40.80
        XNLEN = 0.                                                        40.80
        YNLEN = 0.                                                        40.80
        MXN   = 0                                                         40.80
        MYN   = 0                                                         40.80
        XPCN  = 0.                                                        40.80
        YPCN  = 0.                                                        40.80
        ALPCN = 0.                                                        40.80
      ELSE
        WRITE (PRTEST,'(A)') ' error SWODDC: no PSTYPE defined'           40.31
      ENDIF
!
      IF (PSNAME.EQ.'COMPGRID') THEN
         LCOMPGRD=.TRUE.
      ELSE
         LCOMPGRD=.FALSE.
      ENDIF
!
      IF (ITEST.GE.100 .OR. IOUTES .GE. 30) THEN
        WRITE (PRTEST, 91) PSNAME, PSTYPE, MIP
  91    FORMAT (' Exit SWODDC  ', A16, 2X, A1, 3I6)
        IF (PSTYPE.EQ.'F' .OR. PSTYPE .EQ. 'H') WRITE (PRTEST, 92)
     &     MXK, MYK, ALPQ, DXK, DYK
  92    FORMAT ('SWODDC : MXK, MYK,   ALPQ,   DXK,     DYK',/,
     &  6X, 2I5, 4(1X,E9.3))
      ENDIF
!
      RETURN
!*    end of subroutine SWODDC   **
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOEXC ( PSTYPE    ,OPI        ,OPR   ,                  40.31
     &                    X         ,Y          ,                         40.31
     &                    MIP       ,XP         ,
     &                    YP        ,XC         ,
     &                    YC        ,KGRPNT     ,
     &                    XCGRID    ,YCGRID     ,                         30.21
     &                    CROSS                 )                         40.86 40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde (1D version)
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.86: Nico Booij
!
!  1. Update
!
!     30.72, Sept 97: placed a missing comma in FORMAT statement
!     30.72, Sept 97: Replaced DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     32.02, Feb. 98: Introduced 1D version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, June 98: argument KGRBND added, call CVMESH modified
!     40.02, Oct. 00: Gave KGRBND array dimension
!     40.13, Aug. 01: repeating grid (KREPTX>0) XC modified
!                     swcomm4.inc reactivated
!     40.30, Apr. 03: introduction distributed-memory approach using MPI
!     40.31, Dec. 03: removing POOL mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.86, Feb. 08: modifications to prevent interpolation over obstacles
!                     arguments added to list
!
!  2. Purpose
!
!     Calculates computational grid coordinates of the output points
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
      INTEGER OPI(2)                                                      40.31
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
      REAL    OPR(5), X(MIP), Y(MIP)                                      40.31
!
!     PSTYPE  Char   input    type of output point set
!     MIP     Int    input    number of output points
!     XP, YP  real   outp     user coordinates of output point
!     XC, YC  real   outp     comp. grid coordinates
!
      LOGICAL CROSS(4,MIP) ! true if obstacle is between output point     40.86
                           ! and computational grid point                 40.86
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     OUTPUT (SWAN/OUTP)
!
! 11. Remarks
!
!     ---
!
! 13. Source text
!
      REAL       XC(*), YC(*), XP(MIP), YP(MIP)                           40.31
      CHARACTER  PSTYPE *1
      INTEGER     KGRPNT(MXC,MYC)                                         30.21
      INTEGER   ITMP1, ITMP2, ITMP3, ITMP4, ITMP5, ITMP6
      REAL      RTMP1, RTMP2, RTMP3, RTMP4
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWOEXC')
!
      IF (PSTYPE.EQ.'F') THEN
        MXQ   = OPI(1)                                                    40.31
        MYQ   = OPI(2)                                                    40.31
        ALPQ  = OPR(5)                                                    40.31
        COSPQ = COS(ALPQ)
        SINPQ = SIN(ALPQ)
        XPQ   = OPR(1)                                                    40.31
        YPQ   = OPR(2)                                                    40.31
        XQLEN = OPR(3)                                                    40.31
        YQLEN = OPR(4)                                                    40.31
        XQP   = -XPQ*COSPQ - YPQ*SINPQ
        YQP   =  XPQ*SINPQ - YPQ*COSPQ
        IF (MXQ.GT.1) THEN
          DXQ = XQLEN/(MXQ-1)
        ELSE
          DXQ = 0.01
        ENDIF
        IF (MYQ.GT.1) THEN
          DYQ = YQLEN/(MYQ-1)
        ELSE
          DYQ = 0.01
        ENDIF
        IP    = 0
        DO  11  IYQ = 1, MYQ                                              30.72
          YY  = (IYQ-1)*DYQ
          XP1 = XPQ - YY*SINPQ
          YP1 = YPQ + YY*COSPQ
          DO  10  IXQ = 1, MXQ
            XX = (IXQ-1)*DXQ
            IP = IP+1
            XP(IP) = XP1 + XX*COSPQ
            YP(IP) = YP1 + XX*SINPQ
   10     CONTINUE                                                        30.72
   11   CONTINUE                                                          30.72
      ELSE IF (PSTYPE .EQ. 'H') THEN                                      30.21
        XCMAX = OPR(1)                                                    40.31
        YCMAX = OPR(2)                                                    40.31
        XCMIN = OPR(3)                                                    40.31
        YCMIN = OPR(4)                                                    40.31
        ALPQ  = OPR(5)                                                    40.31
        MXQ   = OPI(1)                                                    40.31
        MYQ   = OPI(2)                                                    40.31
        COSPQ = COS(ALPQ)
        SINPQ = SIN(ALPQ)
        IF (MXQ.GT.1) THEN
           DCXQ = (XCMAX - XCMIN)/(MXQ-1)
        ELSE
           DCXQ = 0.
        END IF
        IF (MYQ.GT.1) THEN
           DCYQ = (YCMAX - YCMIN)/(MYQ-1)
        ELSE
           DCYQ = 0.
        END IF
        XC(1) = XCMIN
        YC(1) = YCMIN
!       *** Find XQLEN and YQLEN taken the extreme points    ***
!       *** that belongs to the frame                        ***
        XPMIN =  1.E09
        YPMIN =  1.E09
        XPMAX = -1.E09
        YPMAX = -1.E09
!
        XPQ   = XPMIN
        YPQ   = YPMIN
        XQP   = 0.
        YQP   = 0.
        XQLEN = XPMAX - XPMIN
        YQLEN = YPMAX - YPMIN
        IF (MXQ.GT.1) THEN
           DXQ = (XQLEN)/(MXQ-1)
        ELSE
           DXQ = 0.
        END IF
        IF (MYQ.GT.1) THEN
           DYQ = (YQLEN)/(MYQ-1)
        ELSE
           DYQ = 0.
        END IF
        IF (ITEST.GE. 120 ) THEN
          WRITE(PRINTF,64) XQLEN ,YQLEN ,MXQ ,MYQ ,
     &                     DCXQ  ,DCYQ ,XC(1) ,YC(1),DXQ ,DYQ
 64       FORMAT (' SWOEXC FRAME DATA :',/,' XQLEN      ,YQLEN      ',    30.72
     &            ',MXQ ,MYQ , DCXQ    ,DCYQ     ,XC(1)    ,YC(1)',
     &            '    ,DXQ       ,DYQ',/,1X,
     &    2(1X,E9.3),2(1X,I4),2X,6(1X,E9.3))
          WRITE(PRINTF,65)XCMIN,XCMAX,YCMIN,YCMAX,
     &                    XPMIN,XPMAX,YPMIN,YPMAX
 65       FORMAT('XCMIN        ,XCMAX  ,YCMIN     ,YCMAX    ,XPMIN',
     &           '   ,XPMAX     ,YPMIN    ,YPMAX   ',/,8(1X,E9.3),/)
        ENDIF
        YY    = YCMIN - DCYQ
        IP    = 0
        DO 15 IYQ = 1 ,MYQ
          XX = XCMIN - DCXQ
          YY = YY    + DCYQ
          DO 16 IXQ = 1 ,MXQ
            IP = IP + 1
            XX = XX + DCXQ
            XC(IP) = XX - REAL(MXF) + 1.                                  40.31
            YC(IP) = YY - REAL(MYF) + 1.                                  40.31
            IF ( XC(IP).GE.-0.01 .AND. XC(IP).LE.REAL(MXC-1)+0.01 .AND.   40.31
     &           YC(IP).GE.-0.01 .AND. YC(IP).LE.REAL(MYC-1)+0.01 ) THEN  40.41 40.31
               IF ( KGRPNT(NINT(XC(IP))+1,NINT(YC(IP))+1).GT.1 ) THEN     40.41 40.31
                  CALL EVALF (XC(IP)+1.,YC(IP)+1.,XPP,YPP,XCGRID,YCGRID)  30.72
                  XP(IP) = XPP
                  YP(IP) = YPP
               ELSE
                  XP(IP) = OVEXCV(1)
                  YP(IP) = OVEXCV(2)
               END IF
            ELSE
              XP(IP) = OVEXCV(1)
              YP(IP) = OVEXCV(2)
            ENDIF
            IF (ITEST.GE.200) WRITE(PRTEST,63)
     &      XP(IP), YP(IP), XC(IP), YC(IP)
 16       CONTINUE
 15     CONTINUE
        GOTO 85                                                           40.86
      ELSE IF (PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P' .OR.                      20.6x
     &         PSTYPE.EQ.'N' .OR. PSTYPE.EQ.'U' ) THEN                    40.80
        XP = X                                                            40.31
        YP = Y                                                            40.31
      ENDIF
!
!     transform to computational grid
!
      IF (ITEST.GE. 150 .AND. OPTG .EQ. 1)
     &  WRITE (PRTEST, 62) XCP, YCP, COSPC, SINPC, DX, DY
  62  FORMAT (' SWOEXC, transf. coeff.:', 8(1X,E12.4))
!     *** The transformation to computational grid depends ***
!     *** on the grid type: regular(1) , curvilinear(3)    ***
      DO 70 IP=1, MIP
        IF (OPTG .EQ. 1) THEN                                             30.2x
          XC(IP) = (XCP + XP(IP)*COSPC + YP(IP)*SINPC) / DX
          XC(IP) = XC(IP) - REAL(MXF) + 1.                                40.30
!         repeating grid: XC is shifted to be between 0 and MXC           40.13
          IF (KREPTX.GT.0) XC(IP) = MODULO (XC(IP), REAL(MXC))            40.13
          IF (ONED) THEN                                                  32.02
            YC(IP) = 0                                                    32.02
          ELSE                                                            32.02
            YC(IP) = (YCP - XP(IP)*SINPC + YP(IP)*COSPC) / DY
            YC(IP) = YC(IP) - REAL(MYF) + 1.                              40.30
          ENDIF                                                           32.02
        ELSEIF (OPTG.EQ.3) THEN                                           40.80
          XPA = XP(IP)
          YPA = YP(IP)
          ITMP1  = MXC
          ITMP2  = MYC
          ITMP3  = MCGRD
          ITMP4  = NGRBND
          ITMP5  = MXF
          ITMP6  = MYF
          RTMP1  = XCLMIN
          RTMP2  = XCLMAX
          RTMP3  = YCLMIN
          RTMP4  = YCLMAX
          MXC    = MXCGL
          MYC    = MYCGL
          MCGRD  = MCGRDGL
          NGRBND = NGRBGL
          MXF    = 1
          MYF    = 1
          XCLMIN = XCGMIN
          XCLMAX = XCGMAX
          YCLMIN = YCGMIN
          YCLMAX = YCGMAX
          CALL CVMESH (XPA, YPA, XCA, YCA, KGRPGL, XGRDGL ,YGRDGL,        30.21
     &                 KGRBGL)                                            40.00
          MXC    = ITMP1
          MYC    = ITMP2
          MCGRD  = ITMP3
          NGRBND = ITMP4
          MXF    = ITMP5
          MYF    = ITMP6
          XCLMIN = RTMP1
          XCLMAX = RTMP2
          YCLMIN = RTMP3
          YCLMAX = RTMP4
          XC(IP) = XCA - REAL(MXF) + 1.                                   40.51
          YC(IP) = YCA - REAL(MYF) + 1.                                   40.51
        ENDIF
        IF (ITEST.GE.250) WRITE(PRTEST,63) XP(IP), YP(IP),
     &                                     XC(IP), YC(IP)
  70  CONTINUE
  63  FORMAT (' SWOEXC, PROBLEM  COORD:', 2(1X,F12.4),/,
     &        '         COMPUT   COORD:', 2(1X,F12.4))
!
!     find crossings of obstacles between output and grid points          40.86
!
  85  CALL SWOBSTO (XCGRID, YCGRID, XP, YP, XC, YC, KGRPNT, CROSS, MIP)   40.86
!
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOEXD (OQPROC, MIP, XC, YC, VOQR, VOQ, COMPDA ,KGRPNT,  30.21
     &                   FORCE, CROSS, IONOD )                            40.86 40.80 40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE TIMECOMM                                                        40.41
      USE M_PARALL                                                        40.31
      USE M_DIFFR                                                         40.21
      USE OUTP_DATA
      USE SwanGriddata                                                    40.80
      USE SwanGridobjects                                                 40.91
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde (1D version)
!     31.02, 40.13: Nico Booij
!     40.21: Agnieszka Herman
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     40.61: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.86: Nico Booij
!     41.12: Nico Booij
!
!  1. Updates
!
!     10.07, July 94, error in current velocity repaired, wind velocity added
!     30.72, Oct. 97: logical function EQREAL introduced for floating point
!                     comparisons
!     32.02, Feb. 98: Introduced 1D version
!     31.02, Sep. 97: computation of Setup, and computation of Force
!                     from setup computation
!     40.00, June 98: Tsec added
!     33.09, Mar. 00: in case of spherical coordinates, distance in m is
!                     calculated from coordinates in degrees
!     40.13, Oct. 01: Forces always computed by SWOEXF
!     40.21, Nov. 01: diffraction parameter added
!     40.41, Aug. 04: friction coefficient added
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: bottom wave period added
!     40.51, Sep. 05: water level and bottom level added
!     40.61, Sep. 06: separate dissipation coefficients added
!     40.80, Sep. 07: extension to unstructured grids
!     40.86, Feb. 08: interpolation near obstacles modified,
!                     points on the other side of the obstacle not taken into account
!                     calls of SWIPOL changed
!     41.12, Apr. 10: output quantity NPL (type nr 70) added
!
!  2. Purpose
!
!     Calculates  Dist, Depth, Ux, Uy, ..
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     FORCE   real   input    wave-induced force                          40.80
!     IONOD   Int    outp     array indicating in which subdomain         40.51
!                             output points are located                   40.51
!     OQPROC  logic  input    y/n process outp quantities
!     PSNAME  Char   input    name of output point set referred to
!     MIP     Int    input    number of output points
!     XP, YP  real   outp     user coordinates of output point
!     XC, YC  real   outp     comp. grid coordinates
!     WX2, WY2  real input    wind components
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTPUT)
!
!  8. Subroutines used
!
!     SWIPOL
!
! 11. Remarks
!
!       ---
!
! 13. Source text
!
      REAL       XC(*), YC(*), VOQ(MIP,*), COMPDA(MCGRD,MCMVAR)
      REAL       FORCE(nverts,2)                                          40.80
      INTEGER    VOQR(*), KGRPNT(MXC,MYC)
      INTEGER    IONOD(*)                                                 40.31
      INTEGER, ALLOCATABLE :: KVERT(:)                                    41.07
      LOGICAL    OQPROC(*), EQREAL                                        30.72
      LOGICAL    CROSS(4,MIP)                                             40.86
      LOGICAL, ALLOCATABLE :: LTMP(:)                                     40.91
      type(verttype), dimension(:), pointer :: vert                       40.91
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWOEXD')
!
      vert => gridobject%vert_grid                                        40.91
      IF (ITEST.GE. 100 .OR. IOUTES .GE. 10) WRITE (PRTEST, 10)
     &(OQPROC(JJ), JJ=1,10), OQPROC(26), MIP
  10  FORMAT (' Entry SWOEXD ', 11L2, I8)
      IVXP   = 1
      IVYP   = 2
!
!     distance
!
 110  IF (OQPROC(3)) THEN
        IVDIST = VOQR(3)
        DO IP = 1, MIP
          IF (IP.EQ.1) THEN
            RDIST = 0.
          ELSE
            RDX = VOQ(IP,IVXP) - VOQ(IP-1,IVXP)
            RDY = VOQ(IP,IVYP) - VOQ(IP-1,IVYP)
            IF (KSPHER.GT.0) THEN                                         33.09
!             spherical coordinates: distance is expressed in m
              RDX = RDX * LENDEG *
     &        COS(DEGRAD*(YOFFS+0.5*(VOQ(IP,IVYP)+VOQ(IP-1,IVYP))))       33.09
              RDY = RDY * LENDEG
            ENDIF
            RDIST = RDIST + SQRT(RDX*RDX+RDY*RDY)
          ENDIF
          VOQ(IP,IVDIST) = RDIST
        ENDDO
      ENDIF
!
      IF (OPTG.EQ.5) THEN
!
!        ---find closest vertex for given point in case of unstructured grid
!
         ALLOCATE(KVERT(MIP))
         IF (.NOT.LCOMPGRD) THEN
            DO IP = 1, MIP
               CALL SwanFindPoint ( VOQ(IP,1), VOQ(IP,2), KVERT(IP) )     41.07
            ENDDO
         ELSE
            KVERT = 0
         ENDIF
!
      ENDIF
!
!     depth
!
 120  IF (OQPROC(4)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 4,
     &  VOQR(4), JDP2
 121    FORMAT (' SWOEXD, type:', 4I3)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JDP2), OVEXCV(4), XC, YC, MIP, CROSS,    40.86
     &                  VOQ(1,VOQR(4)) ,KGRPNT, COMPDA(1,JDP2))
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(4)), VOQ(1,1),         40.80
     &                                  VOQ(1,2), COMPDA(1,JDP2),         40.80
     &                                  MIP, KVERT, OVEXCV(4) )           40.80
        ENDIF                                                             40.80
      ENDIF
!
!     current velocity
!
 130  IF (OQPROC(5)) THEN
        JVQX = VOQR(5)
        JVQY = JVQX+1
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 5,
     &  VOQR(5), JVX2
        IF (ICUR.EQ.1) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL (COMPDA(1,JVX2), OVEXCV(5), XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,JVQX) ,KGRPNT, COMPDA(1,JDP2))            30.21
             CALL SWIPOL (COMPDA(1,JVY2), OVEXCV(5), XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,JVQY) ,KGRPNT, COMPDA(1,JDP2))            30.21
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,JVQX), VOQ(1,1),          40.80
     &                                    VOQ(1,2), COMPDA(1,JVX2),       40.80
     &                                    MIP, KVERT, OVEXCV(5) )         40.80
             CALL SwanInterpolateOutput ( VOQ(1,JVQY), VOQ(1,1),          40.80
     &                                    VOQ(1,2), COMPDA(1,JVY2),       40.80
     &                                    MIP, KVERT, OVEXCV(5) )         40.80
          ENDIF                                                           40.80
          DO IP = 1, MIP
            UXLOC = VOQ(IP,JVQX)
            UYLOC = VOQ(IP,JVQY)
            VOQ(IP,JVQX) = COSCQ*UXLOC - SINCQ*UYLOC
            VOQ(IP,JVQY) = SINCQ*UXLOC + COSCQ*UYLOC
          ENDDO                                                           10.07
        ELSE
          DO IP = 1, MIP                                                  20.85
            VOQ(IP,JVQX) = 0.
            VOQ(IP,JVQY) = 0.
          ENDDO
        ENDIF
      ENDIF
!
!     Ubot
!
 140  IF (OQPROC(6)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 6,
     &  VOQR(6)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JUBOT), OVEXCV(6), XC, YC, MIP, CROSS,   40.86
     &                  VOQ(1,VOQR(6)) ,KGRPNT, COMPDA(1,JDP2))           30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(6)), VOQ(1,1),         40.80
     &                                  VOQ(1,2), COMPDA(1,JUBOT),        40.80
     &                                  MIP, KVERT, OVEXCV(6) )           40.80
        ENDIF                                                             40.80
        KK = VOQR(6)
        RR = SQRT(2.)
        DO IP = 1, MIP
          UBLOC = VOQ(IP,KK)
          IF (.NOT.EQREAL(UBLOC,OVEXCV(6))) VOQ(IP,KK) = RR * UBLOC       30.72
        ENDDO
      ENDIF
!
!     Urms
!
      IF (OQPROC(34)) THEN                                                20.67
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 34,
     &  VOQR(34)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JUBOT), OVEXCV(34), XC, YC, MIP, CROSS,  40.86
     &                  VOQ(1,VOQR(34)) ,KGRPNT, COMPDA(1,JDP2))          30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(34)), VOQ(1,1),        40.80
     &                                  VOQ(1,2), COMPDA(1,JUBOT),        40.80
     &                                  MIP, KVERT, OVEXCV(34) )          40.80
        ENDIF                                                             40.80
      ENDIF
!
!     TmBot
!
      IF (OQPROC(50)) THEN                                                40.51
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 50,
     &  VOQR(50), JPBOT
        IF (JPBOT.GT.1) THEN                                              40.65
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JPBOT),OVEXCV(50),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(50)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(50)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JPBOT),     40.80
     &                                     MIP, KVERT, OVEXCV(50) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(50)) = OVEXCV(50)                                40.65
           ENDDO                                                          40.65
        ENDIF
      ENDIF
!
!     dissipation
!
      IF (OQPROC(7)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 7,
     &  VOQR(7)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JDISS), OVEXCV(7), XC, YC, MIP, CROSS,   40.86
     &                  VOQ(1,VOQR(7)) ,KGRPNT, COMPDA(1,JDP2))           30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(7)), VOQ(1,1),         40.80
     &                                  VOQ(1,2), COMPDA(1,JDISS),        40.80
     &                                  MIP, KVERT, OVEXCV(7) )           40.80
        ENDIF                                                             40.80
        IF (INRHOG.EQ.1) THEN
          DO 152 IP = 1, MIP
            F1 = VOQ(IP,VOQR(7))
            IF (.NOT.EQREAL(F1,OVEXCV(7))) VOQ(IP,VOQR(7))=F1*RHO*GRAV    30.72
 152      CONTINUE
        ENDIF
      ENDIF
!
!     bottom friction dissipation
!
      IF (OQPROC(54)) THEN                                                40.61
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 54,
     &  VOQR(54), JDSXB
        IF (JDSXB.GT.1) THEN                                              40.65
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JDSXB),OVEXCV(54),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(54)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(54)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JDSXB),     40.80
     &                                     MIP, KVERT, OVEXCV(54) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(54)) = OVEXCV(54)                                40.65
           ENDDO                                                          40.65
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(54))
            IF (.NOT.EQREAL(F1,OVEXCV(54))) VOQ(IP,VOQR(54))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     wave breaking dissipation
!
      IF (OQPROC(55)) THEN                                                40.61
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 55,
     &  VOQR(55), JDSXS
        IF (JDSXS.GT.1) THEN                                              40.65
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JDSXS),OVEXCV(55),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(55)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(55)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JDSXS),     40.80
     &                                     MIP, KVERT, OVEXCV(55) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(55)) = OVEXCV(55)                                40.65
           ENDDO                                                          40.65
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(55))
            IF (.NOT.EQREAL(F1,OVEXCV(55))) VOQ(IP,VOQR(55))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     whitecapping dissipation
!
      IF (OQPROC(56)) THEN                                                40.61
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 56,
     &  VOQR(56), JDSXW
        IF (JDSXW.GT.1) THEN                                              40.65
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JDSXW),OVEXCV(56),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(56)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(56)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JDSXW),     40.80
     &                                     MIP, KVERT, OVEXCV(56) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(56)) = OVEXCV(56)                                40.65
           ENDDO                                                          40.65
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(56))
            IF (.NOT.EQREAL(F1,OVEXCV(56))) VOQ(IP,VOQR(56))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     vegetation dissipation
!
      IF (OQPROC(57)) THEN                                                40.61
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 57,
     &  VOQR(57), JDSXV
        IF (JDSXV.GT.1) THEN                                              40.65
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JDSXV),OVEXCV(57),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(57)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(57)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JDSXV),     40.80
     &                                     MIP, KVERT, OVEXCV(57) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(57)) = OVEXCV(57)                                40.65
           ENDDO                                                          40.65
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(57))
            IF (.NOT.EQREAL(F1,OVEXCV(57))) VOQ(IP,VOQR(57))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     energy generation
!
      IF (OQPROC(60)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 60,
     &  VOQR(60), JGENR
        IF (JGENR.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JGENR),OVEXCV(60),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(60)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(60)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JGENR),
     &                                     MIP, KVERT, OVEXCV(60) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(60)) = OVEXCV(60)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(60))
            IF (.NOT.EQREAL(F1,OVEXCV(60))) VOQ(IP,VOQR(60))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     wind source term
!
      IF (OQPROC(61)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 61,
     &  VOQR(61), JGSXW
        IF (JGSXW.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JGSXW),OVEXCV(61),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(61)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(61)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JGSXW),
     &                                     MIP, KVERT, OVEXCV(61) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(61)) = OVEXCV(61)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(61))
            IF (.NOT.EQREAL(F1,OVEXCV(61))) VOQ(IP,VOQR(61))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     energy redistribution
!
      IF (OQPROC(62)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 62,
     &  VOQR(62), JREDS
        IF (JREDS.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JREDS),OVEXCV(62),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(62)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(62)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JREDS),
     &                                     MIP, KVERT, OVEXCV(62) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(62)) = OVEXCV(62)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(62))
            IF (.NOT.EQREAL(F1,OVEXCV(62))) VOQ(IP,VOQR(62))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     total absolute 4-wave interaction
!
      IF (OQPROC(63)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 63,
     &  VOQR(63), JRSXQ
        IF (JRSXQ.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JRSXQ),OVEXCV(63),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(63)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(63)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JRSXQ),
     &                                     MIP, KVERT, OVEXCV(63) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(63)) = OVEXCV(63)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(63))
            IF (.NOT.EQREAL(F1,OVEXCV(63))) VOQ(IP,VOQR(63))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     total absolute 3-wave interaction
!
      IF (OQPROC(64)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 64,
     &  VOQR(64), JRSXT
        IF (JRSXT.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JRSXT),OVEXCV(64),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(64)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(64)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JRSXT),
     &                                     MIP, KVERT, OVEXCV(64) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(64)) = OVEXCV(64)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(64))
            IF (.NOT.EQREAL(F1,OVEXCV(64))) VOQ(IP,VOQR(64))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     energy propagation
!
      IF (OQPROC(65)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 65,
     &  VOQR(65), JTRAN
        IF (JTRAN.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JTRAN),OVEXCV(65),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(65)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(65)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JTRAN),
     &                                     MIP, KVERT, OVEXCV(65) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(65)) = OVEXCV(65)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(65))
            IF (.NOT.EQREAL(F1,OVEXCV(65))) VOQ(IP,VOQR(65))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     xy-propagation
!
      IF (OQPROC(66)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 66,
     &  VOQR(66), JTSXG
        IF (JTSXG.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JTSXG),OVEXCV(66),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(66)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(66)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JTSXG),
     &                                     MIP, KVERT, OVEXCV(66) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(66)) = OVEXCV(66)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(66))
            IF (.NOT.EQREAL(F1,OVEXCV(66))) VOQ(IP,VOQR(66))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     theta-propagation
!
      IF (OQPROC(67)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 67,
     &  VOQR(67), JTSXT
        IF (JTSXT.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JTSXT),OVEXCV(67),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(67)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(67)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JTSXT),
     &                                     MIP, KVERT, OVEXCV(67) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(67)) = OVEXCV(67)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(67))
            IF (.NOT.EQREAL(F1,OVEXCV(67))) VOQ(IP,VOQR(67))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     sigma-propagation
!
      IF (OQPROC(68)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 68,
     &  VOQR(68), JTSXS
        IF (JTSXS.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JTSXS),OVEXCV(68),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(68)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(68)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JTSXS),
     &                                     MIP, KVERT, OVEXCV(68) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(68)) = OVEXCV(68)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(68))
            IF (.NOT.EQREAL(F1,OVEXCV(68))) VOQ(IP,VOQR(68))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     radiation stress
!
      IF (OQPROC(69)) THEN                                                40.85
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 69,
     &  VOQR(69), JRADS
        IF (JRADS.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JRADS),OVEXCV(69),XC, YC, MIP, CROSS,
     &                    VOQ(1,VOQR(69)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(69)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JRADS),
     &                                     MIP, KVERT, OVEXCV(69) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(69)) = OVEXCV(69)
           ENDDO
        ENDIF
        IF (INRHOG.EQ.1) THEN
          DO IP = 1, MIP
            F1 = VOQ(IP,VOQR(69))
            IF (.NOT.EQREAL(F1,OVEXCV(69))) VOQ(IP,VOQR(69))=F1*RHO*GRAV
          END DO
        ENDIF
      ENDIF
!
!     number of plants per square meter
!
      IF (OQPROC(70)) THEN                                              ! 41.12
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 70,
     &  VOQR(70), JNPLA2
        IF (JNPLA2.GT.1) THEN
           IF (OPTG.NE.5) THEN
              CALL SWIPOL(COMPDA(1,JNPLA2),OVEXCV(70),XC,YC, MIP, CROSS,
     &                    VOQ(1,VOQR(70)) ,KGRPNT, COMPDA(1,JDP2))
           ELSE
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(70)), VOQ(1,1),
     &                                     VOQ(1,2), COMPDA(1,JNPLA2),
     &                                     MIP, KVERT, OVEXCV(70) )
           ENDIF
        ELSE
           DO IP = 1, MIP
             VOQ(IP,VOQR(70)) = OVEXCV(70)
           ENDDO
        ENDIF
      ENDIF
!
!     Qb
!
      IF (OQPROC(8)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 8,
     &  VOQR(8)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JQB), OVEXCV(8), XC, YC, MIP, CROSS,     40.86
     &                  VOQ(1,VOQR(8)) ,KGRPNT, COMPDA(1,JDP2))           30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(8)), VOQ(1,1),         40.80
     &                                  VOQ(1,2), COMPDA(1,JQB),          40.80
     &                                  MIP, KVERT, OVEXCV(8) )           40.80
        ENDIF                                                             40.80
      ENDIF
!
!     wind velocity                                                       10.07
!
      IF (OQPROC(26)) THEN
        JVQX = VOQR(26)
        JVQY = JVQX+1
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 26,
     &  VOQR(26)
        IF (VARWI) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL (COMPDA(1,JWX2), OVEXCV(26),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,JVQX) ,KGRPNT, COMPDA(1,JDP2))
             CALL SWIPOL (COMPDA(1,JWY2), OVEXCV(26),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,JVQY) ,KGRPNT, COMPDA(1,JDP2))            30.21
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,JVQX), VOQ(1,1),          40.80
     &                                    VOQ(1,2), COMPDA(1,JWX2),       40.80
     &                                    MIP, KVERT, OVEXCV(26) )        40.80
             CALL SwanInterpolateOutput ( VOQ(1,JVQY), VOQ(1,1),          40.80
     &                                    VOQ(1,2), COMPDA(1,JWY2),       40.80
     &                                    MIP, KVERT, OVEXCV(26) )        40.80
          ENDIF                                                           40.80
          DO IP = 1, MIP
            UXLOC = VOQ(IP,JVQX)
            UYLOC = VOQ(IP,JVQY)
            VOQ(IP,JVQX) = COSCQ*UXLOC - SINCQ*UYLOC
            VOQ(IP,JVQY) = SINCQ*UXLOC + COSCQ*UYLOC
          ENDDO
        ELSE
          DO IP = 1, MIP
            VOQ(IP,JVQX) = U10*COS(WDIP-ALPQ)                             10.36
            VOQ(IP,JVQY) = U10*SIN(WDIP-ALPQ)                             10.36
          ENDDO
        ENDIF
      ENDIF
!
!     difference in Hs between iterations                                 20.52
!
 180  IF (OQPROC(30)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 30,
     &  VOQR(30), JDHS
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JDHS), OVEXCV(30), XC, YC, MIP, CROSS,   40.86
     &                  VOQ(1,VOQR(30)) ,KGRPNT, COMPDA(1,JDP2))          30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(30)), VOQ(1,1),        40.80
     &                                  VOQ(1,2), COMPDA(1,JDHS),         40.80
     &                                  MIP, KVERT, OVEXCV(30) )          40.80
        ENDIF                                                             40.80
      ENDIF
!
!     difference in Tm between iterations                                 20.52
!
 190  IF (OQPROC(31)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 31,
     &  VOQR(31), JDTM
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JDTM), OVEXCV(31), XC, YC, MIP, CROSS,   40.86
     &                  VOQ(1,VOQR(31)) ,KGRPNT, COMPDA(1,JDP2))          30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(31)), VOQ(1,1),        40.80
     &                                  VOQ(1,2), COMPDA(1,JDTM),         40.80
     &                                  MIP, KVERT, OVEXCV(31) )          40.80
        ENDIF                                                             40.80
      ENDIF
!
!     leak
!
 200  IF (OQPROC(9)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 9,
     &  VOQR(9), JLEAK
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JLEAK), OVEXCV(9), XC, YC, MIP, CROSS,   40.86
     &                  VOQ(1,VOQR(9)) ,KGRPNT, COMPDA(1,JDP2))           30.21
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(9)), VOQ(1,1),         40.80
     &                                  VOQ(1,2), COMPDA(1,JLEAK),        40.80
     &                                  MIP, KVERT, OVEXCV(9) )           40.80
        ENDIF                                                             40.80
        IF (INRHOG.EQ.1) THEN
          DO 202 IP = 1, MIP
            F1 = VOQ(IP,VOQR(9))
            IF (.NOT.EQREAL(F1,OVEXCV(9))) VOQ(IP,VOQR(9))=F1*RHO*GRAV    30.72
 202      CONTINUE
        ENDIF
      ENDIF
!
!     Ufric
!
      IF (OQPROC(35)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 35,
     &  VOQR(35), JUSTAR
        IF (JUSTAR.GT.1) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL(COMPDA(1,JUSTAR),OVEXCV(35),XC, YC, MIP, CROSS,  40.86
     &                   VOQ(1,VOQR(35)) ,KGRPNT, COMPDA(1,JDP2))         30.22
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(35)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), COMPDA(1,JUSTAR),     40.80
     &                                    MIP, KVERT, OVEXCV(35) )        40.80
          ENDIF                                                           40.80
        ELSE
          DO IP = 1, MIP                                                  31.02
            VOQ(IP,VOQR(35)) = OVEXCV(35)                                 31.02
          ENDDO                                                           31.02
        ENDIF
      ENDIF
!
!     zelen
!
      IF (OQPROC(36)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 36,
     &  VOQR(36), JZEL
        IF (JZEL.GT.1) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL(COMPDA(1,JZEL), OVEXCV(36), XC, YC, MIP, CROSS,  40.86
     &                   VOQ(1,VOQR(36)) ,KGRPNT, COMPDA(1,JDP2))         30.22
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(36)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), COMPDA(1,JZEL),       40.80
     &                                    MIP, KVERT, OVEXCV(36) )        40.80
          ENDIF                                                           40.80
        ELSE
          DO IP = 1, MIP                                                  31.02
            VOQ(IP,VOQR(36)) = OVEXCV(36)                                 31.02
          ENDDO                                                           31.02
        ENDIF
      ENDIF
!
!     TauW
!
      IF (OQPROC(37)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 37,
     &  VOQR(37), JTAUW
        IF (JTAUW.GT.1) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL(COMPDA(1,JTAUW),OVEXCV(37), XC, YC, MIP, CROSS,  40.86
     &                   VOQ(1,VOQR(37)) ,KGRPNT, COMPDA(1,JDP2))         30.22
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(37)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), COMPDA(1,JTAUW),      40.80
     &                                    MIP, KVERT, OVEXCV(37) )        40.80
          ENDIF                                                           40.80
        ELSE
          DO IP = 1, MIP                                                  31.02
            VOQ(IP,VOQR(37)) = OVEXCV(37)                                 31.02
          ENDDO                                                           31.02
        ENDIF
      ENDIF
!
!     Cdrag
!
      IF (OQPROC(38)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 20) WRITE (PRTEST, 121) 38,
     &  VOQR(38), JCDRAG
        IF (JCDRAG.GT.1) THEN
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL(COMPDA(1,JCDRAG),OVEXCV(38),XC, YC, MIP, CROSS,  40.86
     &                   VOQ(1,VOQR(38)) ,KGRPNT, COMPDA(1,JDP2))         30.22
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(38)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), COMPDA(1,JCDRAG),     40.80
     &                                    MIP, KVERT, OVEXCV(38) )        40.80
          ENDIF                                                           40.80
        ELSE
          DO IP = 1, MIP                                                  31.02
            VOQ(IP,VOQR(38)) = OVEXCV(38)                                 31.02
          ENDDO                                                           31.02
        ENDIF
      ENDIF
!
!     wave-induced setup                                                  32.02
!
      IF (OQPROC(39)) THEN                                                32.02
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 39,      32.02
     &  VOQR(39), JSETUP                                                  32.02
        IF (LSETUP.GT.0) THEN                                             32.02
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL(COMPDA(1,JSETUP),OVEXCV(39),XC, YC, MIP, CROSS,  40.86 32.02
     &                   VOQ(1,VOQR(39)) ,KGRPNT, COMPDA(1,JDP2))         32.02
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(39)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), COMPDA(1,JSETUP),     40.80
     &                                    MIP, KVERT, OVEXCV(39) )        40.80
          ENDIF                                                           40.80
        ELSE                                                              32.02
          DO IP = 1, MIP                                                  32.02
            VOQ(IP,VOQR(39)) = OVEXCV(39)                                 32.02
          ENDDO                                                           32.02
        ENDIF                                                             32.02
      ENDIF                                                               32.02
!
!     wave-induced force (unstructured grids only!)                       40.80
!
      IF (OQPROC(20).AND.OPTG.EQ.5) THEN                                  40.80
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 20,      40.80
     &  VOQR(20)                                                          40.80
        CALL SwanInterpolateOutput ( VOQ(1,VOQR(20)), VOQ(1,1),           40.80
     &                               VOQ(1,2), FORCE(1,1),                40.80
     &                               MIP, KVERT, OVEXCV(20) )             40.80
        CALL SwanInterpolateOutput ( VOQ(1,VOQR(20)+1), VOQ(1,1),         40.80
     &                               VOQ(1,2), FORCE(1,2),                40.80
     &                               MIP, KVERT, OVEXCV(20) )             40.80
      ENDIF                                                               40.80
!
!     Ursell
!
      IF (OQPROC(45)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 45,      40.03
     &  VOQR(45)
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWIPOL (COMPDA(1,JURSEL), OVEXCV(45),XC, YC, MIP, CROSS,  40.86
     &                  VOQ(1,VOQR(45)) ,KGRPNT, COMPDA(1,JDP2))          40.03
        ELSE                                                              40.80
           CALL SwanInterpolateOutput ( VOQ(1,VOQR(45)), VOQ(1,1),        40.80
     &                                  VOQ(1,2), COMPDA(1,JURSEL),       40.80
     &                                  MIP, KVERT, OVEXCV(45) )          40.80
        ENDIF                                                             40.80
      ENDIF
!
!     Air-Sea temperature difference
!
      IF (OQPROC(46)) THEN
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 46,      40.03
     &  VOQR(46), JASTD2
        IF (VARAST) THEN
           IF (OPTG.NE.5) THEN                                            40.80
              CALL SWIPOL(COMPDA(1,JASTD2),OVEXCV(46),XC,YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(46)) ,KGRPNT, COMPDA(1,JDP2))        40.03
           ELSE                                                           40.80
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(46)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JASTD2),    40.80
     &                                     MIP, KVERT, OVEXCV(46) )       40.80
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP
              VOQ(IP,VOQR(46)) = OVEXCV(46)
           END DO
        END IF
      ENDIF
!
!       Diffraction parameter                                             40.21
!
      IF (OQPROC(49)) THEN
        IF (IDIFFR.EQ.1) THEN
          IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 49,    40.21
     &    VOQR(49)
          IF (OPTG.NE.5) THEN                                             40.80
             CALL SWIPOL (DIFPARAM(:), OVEXCV(49), XC, YC, MIP, CROSS,    40.86 40.21
     &                    VOQ(1,VOQR(49)) ,KGRPNT, COMPDA(1,JDP2))        40.21
          ELSE                                                            40.80
             CALL SwanInterpolateOutput ( VOQ(1,VOQR(49)), VOQ(1,1),      40.80
     &                                    VOQ(1,2), DIFPARAM(:),          40.80
     &                                    MIP, KVERT, OVEXCV(49) )        40.80
          ENDIF                                                           40.80
        ELSE
          DO IP = 1, MIP                                                  40.21
            VOQ(IP,VOQR(49)) = 1.                                         40.21
          ENDDO                                                           40.21
        ENDIF
      ENDIF
!
!     Tsec
!
      IF (OQPROC(41)) THEN
        DO IP = 1, MIP
          VOQ(IP,VOQR(41)) = TIMCO - OUTPAR(1)                            40.00
        ENDDO
      ENDIF
!
!     friction coefficient                                                40.41
!
      IF (OQPROC(27)) THEN
         IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 27,
     &   VOQR(27)
         IF (VARFR) THEN
            IF (OPTG.NE.5) THEN                                           40.80
!              interpolation done in all active and non-active points
               RTMP=DEPMIN
               DEPMIN=-10.*ABS(MINVAL(COMPDA(:,JDP2)))
               CALL SWIPOL(COMPDA(1,JFRC2),OVEXCV(27),XC,YC, MIP, CROSS,  40.86
     &                     VOQ(1,VOQR(27)) ,KGRPNT, COMPDA(1,JDP2))
               DEPMIN=RTMP
            ELSE                                                          40.80
!              interpolation done in all active and non-active points     40.91
               ALLOCATE(LTMP(nverts))                                     40.91
               LTMP(:) = vert(:)%active                                   40.91
               vert(:)%active = .TRUE.                                    40.91
               CALL SwanInterpolateOutput ( VOQ(1,VOQR(27)), VOQ(1,1),    40.80
     &                                      VOQ(1,2), COMPDA(1,JFRC2),    40.80
     &                                      MIP, KVERT, OVEXCV(27) )      40.80
               vert(:)%active = LTMP(:)                                   40.91
               DEALLOCATE(LTMP)                                           40.91
            ENDIF                                                         40.80
         ELSE
            F1=0.
            IF (IBOT.EQ.1) F1 = PBOT(3)
            IF (IBOT.EQ.2) F1 = PBOT(2)
            IF (IBOT.EQ.3) F1 = PBOT(5)
            DO IP = 1, MIP
               IF (.NOT.EQREAL(F1,OVEXCV(27))) VOQ(IP,VOQR(27))=F1
            END DO
        END IF
      ENDIF
!
!     water level
!
      IF (OQPROC(51)) THEN                                                40.51
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 51,
     &  VOQR(51)
        IF (VARWLV) THEN
           IF (OPTG.NE.5) THEN                                            40.80
!             interpolation done in all active and non-active points
              RTMP=DEPMIN
              DEPMIN=-10.*ABS(MINVAL(COMPDA(:,JDP2)))
              CALL SWIPOL(COMPDA(1,JWLV2),OVEXCV(51),XC, YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(51)) ,KGRPNT, COMPDA(1,JDP2))
              DEPMIN=RTMP
           ELSE                                                           40.80
!             interpolation done in all active and non-active points      40.91
              ALLOCATE(LTMP(nverts))                                      40.91
              LTMP(:) = vert(:)%active                                    40.91
              vert(:)%active = .TRUE.                                     40.91
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(51)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JWLV2),     40.80
     &                                     MIP, KVERT, OVEXCV(51) )       40.80
              vert(:)%active = LTMP(:)                                    40.91
              DEALLOCATE(LTMP)                                            40.91
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP
              VOQ(IP,VOQR(51)) = 0.
           END DO
        END IF
        DO IP = 1, MIP
           F1 = VOQ(IP,VOQR(51))
           IF (.NOT.EQREAL(F1,OVEXCV(51))) VOQ(IP,VOQR(51))=F1 + WLEV
        END DO
      ENDIF
!
!     bottom level
!
      IF (OQPROC(52)) THEN                                                40.51
        IF (ITEST.GE.50 .OR. IOUTES .GE. 10) WRITE (PRTEST, 121) 52,
     &  VOQR(52), JBOTLV
        IF (JBOTLV.GT.1) THEN                                             40.65
           IF (OPTG.NE.5) THEN                                            40.80
!             interpolation done in all active and non-active points
              RTMP=DEPMIN
              DEPMIN=-10.*ABS(MINVAL(COMPDA(:,JDP2)))
              CALL SWIPOL(COMPDA(1,JBOTLV),OVEXCV(52),XC,YC, MIP, CROSS,  40.86
     &                    VOQ(1,VOQR(52)) ,KGRPNT, COMPDA(1,JDP2))
              DEPMIN=RTMP
           ELSE                                                           40.80
!             interpolation done in all active and non-active points      40.91
              ALLOCATE(LTMP(nverts))                                      40.91
              LTMP(:) = vert(:)%active                                    40.91
              vert(:)%active = .TRUE.                                     40.91
              CALL SwanInterpolateOutput ( VOQ(1,VOQR(52)), VOQ(1,1),     40.80
     &                                     VOQ(1,2), COMPDA(1,JBOTLV),    40.80
     &                                     MIP, KVERT, OVEXCV(52) )       40.80
              vert(:)%active = LTMP(:)                                    40.91
              DEALLOCATE(LTMP)                                            40.91
           ENDIF                                                          40.80
        ELSE
           DO IP = 1, MIP                                                 40.65
             VOQ(IP,VOQR(52)) = OVEXCV(52)                                40.65
           ENDDO                                                          40.65
        ENDIF
      ENDIF
!
!     correct problem coordinates with offset values
!
      DO IP=1, MIP
        XP1 = VOQ(IP,IVXP)
        IF (.NOT.EQREAL(XP1,OVEXCV(1)))                                   40.00
     &        VOQ(IP,IVXP) = XP1 + XOFFS
        YP1 = VOQ(IP,IVYP)
        IF (.NOT.EQREAL(YP1,OVEXCV(2)))                                   40.00
     &        VOQ(IP,IVYP) = YP1 + YOFFS
      ENDDO
!
!     --- in case of parallel run, mark location points inside own        40.31
!         subdomain                                                       40.31

      IF ( PARLL ) THEN                                                   40.31
         IXB = 1+IHALOX                                                   40.31
         IF ( LMXF ) IXB = 1                                              40.41 40.31
         IXE = MXC-IHALOX                                                 40.31
         IF ( LMXL ) IXE = MXC                                            40.41 40.31
         IYB = 1+IHALOY                                                   40.31
         IF ( LMYF ) IYB = 1                                              40.41 40.31
         IYE = MYC-IHALOY                                                 40.31
         IF ( LMYL ) IYE = MYC                                            40.41 40.31
         DO IP = 1, MIP                                                   40.31
            IX = NINT(XC(IP)+100.) - 99                                   41.07 40.31
            IY = NINT(YC(IP)+100.) - 99                                   41.07 40.31
            IF ( IX.GE.IXB .AND. IX.LE.IXE .AND.                          40.31
     &           IY.GE.IYB .AND. IY.LE.IYE ) IONOD(IP) = INODE            40.31
         END DO                                                           40.31
      END IF                                                              40.31
!
      IF (ALLOCATED(KVERT)) DEALLOCATE(KVERT)                             41.07
!
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWIPOL (FINP, EXCVAL, XC, YC, MIP, CROSS, FOUTP,         40.86
     &                   KGRPNT, DEP2)                                    40.86 40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!     40.86: Nico Booij
!
!  1. UPDATE
!
!     40.00, July 98: no interpolation if one or more corners are dry
!                     argument DEP2 added
!                     margin around comp. grid introduced
!     40.13, Aug. 01: provision for repeating grid
!                     swcomm4.inc reactivated
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.86, Feb. 08: interpolation over an obstacle prevented
!
!  2. PURPOSE
!
!       Interpolate the function FINP to the point given by computational
!       grid coordinates XC and YC; result appears in array FOUTP
!
!  3. METHOD
!
!       This subroutine computes the contributions from surrounding
!       points to the function value in an output point. The points used
!       are indicated in the sketch below.
!
!                                                            Y
!           +-------------------------------------------------->
!           |
!           |       .       .       .       .       .       .
!           |
!           |
!           |
!           |       .       .       *       *       .       .
!           |
!           |
!           |                         o
!           |       .       .       *       *       .       .
!           |
!           |
!           |
!           |       .       .       .       .       .       .
!           |
!           |
!           |
!         X |       .       .       .       .       .       .
!           |
!           V
!
!                 *   point of the computational grid contributing
!                       to output point (o)
!                 .   other grid points
!
!  4. PARAMETERLIST
!
!       FINP    real a input    array of function values defined on the
!                               computational grid
!       EXCVAL  real   input    exception value (assigned if point is outside
!                               computational grid)
!       XC, YC  real a input    array containing computational grid coordinates
!                               of output points
!       MIP     INT    input    number of output points
!       FOUTP   real a output   array of interpolated values for the output
!                               points
!
!  5. SUBROUTINES CALLING
!
!       SWOEXD
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       IINTPC=1: bilinear interpolation
!       IINTPC=2: higher order interpolation using functions G1 and G2
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       For every output point do
!           If the output point is near line XCL then
!               Determine points contributing to the output point
!               Compute contribution to the projection of the output
!                 point on line XCL
!               Compute multiplication factor for interpolation in X-
!                 direction
!               Compute contribution for the output point
!               Add result to value of variable for the output point in
!                 array IFOP
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      REAL FINP(MCGRD), FOUTP(MIP), XC(MIP), YC(MIP), DEP2(MCGRD)         40.00
      LOGICAL CROSS(4,MIP) ! true if obstacle is between output point     40.86
                           ! and computational grid point                 40.86
      LOGICAL OUTSID
      INTEGER  KGRPNT(MXC,MYC)                                            30.21
!
      REAL :: WW(1:4)    ! Interpolation weights for the 4 corners        40.86
      REAL :: SUMWW      ! sum of the weights                             40.86
      INTEGER :: JX(1:4), JY(1:4) ! grid counters for the 4 corners       40.86
      INTEGER :: INDX(1:4)     ! grid counters for the 4 corners          40.86
      INTEGER :: JC            ! corner counter                           40.86
!
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL  STRACE (IENT, 'SWIPOL')
!
        IF (ITEST.GE.150) WRITE (PRTEST, 61)                             060997
  61    FORMAT ('   XC    , YC  ,',
     &  '   JX1, JY1, JX2,  JY2  SX1,  SY1, FOUTP(IP),',
     &  '    INDX1  INDX2  INDX3  INDX4')
!
      DO 100 IP=1,MIP
        IF (XC(IP) .LE. -0.5 .OR. YC(IP) .LE. -0.5) THEN                  40.00
          FOUTP(IP) = EXCVAL
          JX1   = 0
          JY1   = 0
          JX2   = 0
          JY2   = 0
          SX1   = 0.
          SX2   = 0.
          INDX(1:4) = 0                                                   40.86
          GOTO 80
        ENDIF                                                             30.21
        OUTSID = .FALSE.
        FOUTP(IP) = 0.
        JX1 = INT(XC(IP)+3.001) - 2
        JX2 = JX1 + 1
        SX2 = XC(IP) + 1. - FLOAT(JX1)
        SX1 = 1. - SX2
        IF (JX1.LT.0)   OUTSID = .TRUE.
        IF (KREPTX .EQ. 0) THEN                                           40.13
          IF (JX1.GT.MXC) OUTSID = .TRUE.
          IF (JX1.EQ.MXC) JX2 = MXC
          IF (JX1.EQ.0)   JX1 = 1
        ELSE                                                              40.13
          JX1 = 1 + MODULO (JX1-1,MXC)                                    40.13
          JX2 = 1 + MODULO (JX2-1,MXC)                                    40.13
        ENDIF                                                             40.13
        IF (ONED) THEN
          JY1 = 1                                                         40.86
          JY2 = 1                                                         40.86
          SY1 = 0.5                                                       40.86
          SY2 = 0.5                                                       40.86
        ELSE
          JY1 = INT(YC(IP)+3.001) - 2
          JY2 = JY1 + 1
          SY2 = YC(IP) + 1. - FLOAT(JY1)
          SY1 = 1. - SY2
          IF (JY1.LT.0)   OUTSID = .TRUE.
          IF (JY1.GT.MYC) OUTSID = .TRUE.
          IF (JY1.EQ.MYC) JY2 = MYC
          IF (JY1.EQ.0)   JY1 = 1
        ENDIF
        IF (OUTSID) THEN
          FOUTP(IP) = EXCVAL
        ELSE
          JX(1) = JX1                                                     40.86 30.21
          JY(1) = JY1                                                     40.86 30.21
          WW(1) = SX1*SY1                                                 40.86
          JX(2) = JX2                                                     40.86 30.21
          JY(2) = JY1                                                     40.86 30.21
          WW(2) = SX2*SY1                                                 40.86
          JX(3) = JX1                                                     40.86 30.21
          JY(3) = JY2                                                     40.86 30.21
          WW(3) = SX1*SY2                                                 40.86
          JX(4) = JX2                                                     40.86 30.21
          JY(4) = JY2                                                     40.86 30.21
          WW(4) = SX2*SY2                                                 40.86
          DO JC = 1, 4                                                    40.86
            INDX(JC) = KGRPNT(JX(JC),JY(JC))                              40.86 30.21
            IF (WW(JC).LT.0.01) THEN
              WW(JC) = 0.                                                 40.86
            ELSE
              IF (INDX(JC).LE.1) THEN                                     40.86
                WW(JC) = 0.                                               40.86
              ELSE IF (DEP2(INDX(JC)).LE.DEPMIN) THEN                     40.86
                OUTSID = .TRUE.                                           40.94 40.86
              ELSE IF (CROSS(JC,IP) .AND. WW(JC).LT.0.999) THEN           40.86
                WW(JC) = 0.                                               40.86
              ENDIF
            ENDIF
          ENDDO
          SUMWW = SUM(WW(1:4))                                            40.86
          IF (OUTSID) THEN
            FOUTP(IP) = EXCVAL
          ELSE
            IF (SUMWW.GT.0.1) THEN                                        40.86
              FOUTP(IP) = SUM(WW*FINP(INDX)) / SUMWW                      40.86
            ELSE
              FOUTP(IP) = EXCVAL                                          40.86
            ENDIF
          ENDIF
        ENDIF
  80    IF (ITEST.GE.150) WRITE (PRTEST, 82)
     &  XC(IP) , YC(IP) ,JX1, JY1, JX2,JY2, (WW(JC),JC=1,4),              40.86
     &  (INDX(JC), JC=1,4), (FINP(INDX(JC)), JC=1,4)                      40.86
  82    FORMAT (2(F7.1,1X),4I5, 4(1X,F5.2), 3X,4(2X,I5), 3X,              40.86
     &  4(1X,E9.3))                                                       40.86
 100  CONTINUE
!
      RETURN
! * end of subroutine SWIPOL *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOEXA (OQPROC     ,BKC        ,
     &                   MIP        ,XC         ,
     &                   YC         ,VOQR       ,
     &                   VOQ        ,AC2        ,
     &                   ACLOC      ,SPCSIG     ,                         30.72
     &                   WK         ,CG         ,
     &                   SPCDIR     ,NE         ,
     &                   NED        ,KGRPNT     ,
     &                   DEPXY      ,CROSS      )                         40.86 30.50
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.80
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.70, 40.13: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     32.01: Roeland Ris & Cor van der Schelde
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.86: Nico Booij
!     40.87: Marcel Zijlema
!
!  1. Updates
!
!     10.09, Aug. 94: relative and absolute period distinguished
!                     NVOTP increased and type 28 added
!     10.10, Aug. 94: arrays ECOS, ESIN, NE and NED added to arg. list
!     10.22, Sep. 94: condition for tail changed from MSC.GE.3 to MSC.GT.3
!     20.59, Sep. 95: average wave number can be determined with
!                     other powers of k (i.e. OUTPAR(3))
!     20.61, Sep. 95: Tm02 and FWID added; computation of average period
!                     also changed
!     30.72, Oct. 97: logical function EQREAL introduced for floating point
!                     comparisons
!     32.01, Jan. 98: Nautical convention introduced (project h3268)
!     30.70, Feb. 98: ALCQ ignored if nautical direction is requested
!                     computation of kappa corrected (power of Sigma)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Dec. 98: Argument list KSCIP1 adjusted
!     40.13, Aug. 01: provision for repeating grid (KREPTX>0)
!     40.30, May  03: introduction distributed-memory approach using MPI
!     40.41, Sep. 04: added Tm-10 and RTm-10
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 06: added Tp based on parabolic fitting
!     40.80, Sep. 07: extension to unstructured grids
!     40.86, Feb. 08: modification to prevent interpolation over an obstacle
!     40.87, Apr. 08: integration over [fmin,fmax] added
!
!  2. Purpose
!
!     calculates quantities for which the spectral action density is
!     necessary
!
!  3. Method
!
!       ---
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
!
!     OQPROC  logic  input    processing of output quantities
!     MIP     Int    input    number of output points
!     XC, YC  real   input    comp. grid coordinates
!     VOQR    int a  input    location in VOQ of a certain outp quant.
!     VOQ     real a outp     values of output quantities
!     AC2     real a input    action densities
!     WK      real a local    wavenumber in output point
!     CG      real a local    group velocity in output point
!     NE      real a local    ratio of group and phase velocity
!     NED     real a local    derivative of NE with respect to depth
!     DEPXY   real a input    depth in points of computational grid
!
!     Local Variables
!
!     IVOTP   type indicators of output quantities processed by this subr.
!             used for assignment of exception values
!
!  8. Subroutines used
!
!     DEGCNV: Transforms dir. from nautical to cartesian or vice versa    32.01
!     ANGDEG: Transforms degrees to radians                               32.01
!     SWOINA: interpolates 2D action density spectrum
!
!  9. Subroutines calling
!
!     OUTPUT (SWAN/OUTP)
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     For all output points do
!         interpolate action density to the output point
!         If processing of the quantity is requested
!         Then compute wave height, period etc.
!         ------------------------------------------------------------
!         If computation of Cg and K is necessary
!         Then get depth from array VOQ
!              Call KSCIP1 (computes Cg and K)
!         ------------------------------------------------------------
!         If processing of the quantity is requested
!         Then compute energy transport, wavelength etc.
!     ----------------------------------------------------------------
!
! 13. Source text
!
      PARAMETER  (NVOTP=22)                                               41.15 40.64 40.51 40.41 40.00
      REAL       XC(MIP)        ,YC(MIP)       ,AC2(MDC,MSC,MCGRD),
     &           VOQ(MIP,*)     ,
     &           WK(*)          ,
     &           CG(*)          ,ACLOC(MDC,MSC),
     &           NE(*)          ,NED(*)        ,DEPXY(MCGRD), ECS(MDC)
      LOGICAL    CROSS(4,MIP)                                             40.86
!
      INTEGER    VOQR(*)        ,BKC           ,IVOTP(NVOTP)      ,
     &           KGRPNT(MXC,MYC)                                          30.21
!
      LOGICAL    OQPROC(*), EQREAL                                        30.72
      LOGICAL :: EXCPT     ! if true value in point is undefined          40.86
      SAVE IENT, IVOTP
      DATA IENT /0/
      DATA IVOTP /10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 28, 32, 33,     20.61
     &            42, 43, 44, 47, 48, 53, 58, 59, 71/                     41.15 40.64 40.51 40.41 40.00
      CALL STRACE (IENT, 'SWOEXA')
!
!     loop over all output points
!
      DO 800 IP=1,MIP
        DEP = VOQ(IP,VOQR(4))
!
!       assign exception value if depth is negative or point is outside grid
!
        IF (DEP.LE.0.)                    GOTO 700
        IF (EQREAL(DEP,OVEXCV(4)))        GOTO 700                        30.72
        IF (OPTG.NE.5) THEN                                               40.80
           IF (KREPTX.EQ.0) THEN                                          40.13
!            non-repeating grid                                           40.13
             IF (XC(IP) .LT. -0.01)            GOTO 700
             IF (XC(IP) .GT. REAL(MXC-1)+0.01) GOTO 700
           ENDIF                                                          40.13
           IF (YC(IP) .LT. -0.01)            GOTO 700
           IF (YC(IP) .GT. REAL(MYC-1)+0.01) GOTO 700
        ENDIF                                                             40.80
!
!       first the action density spectrum is interpolated
!
        IF (OPTG.NE.5) THEN                                               40.80
           CALL SWOINA (XC(IP), YC(IP), AC2, ACLOC, KGRPNT, DEPXY,        40.86 30.50
     &                  CROSS(1,IP), EXCPT)                               40.86
        ELSE                                                              40.80
           IF (.NOT.LCOMPGRD) THEN
              IF (.NOT.EQREAL(VOQ(IP,1),OVEXCV(1))) XP=VOQ(IP,1)-XOFFS    40.80
              IF (.NOT.EQREAL(VOQ(IP,2),OVEXCV(2))) YP=VOQ(IP,2)-YOFFS    40.80
              CALL SwanInterpolateAc ( ACLOC, XP, YP, AC2, EXCPT )        40.80
           ELSE
              ACLOC(:,:) = AC2(:,:,IP)
           ENDIF
        ENDIF                                                             40.80
        IF (EXCPT) GOTO 700                                               40.86
!
!       Coefficients for high frequency tail
!
        EFTAIL = 1. / (PWTAIL(1) - 1.)
!
!       significant wave height
!
        IVTYPE = 10
        IF (OQPROC(IVTYPE)) THEN
          IF (OUTPAR(6).EQ.0.) THEN                                       40.87
!            integration over [0,inf]                                     40.87
             ETOT = 0.
!            trapezoidal rule is applied
             DO ID=1, MDC
               DO IS=2,MSC
                 DS=SPCSIG(IS)-SPCSIG(IS-1)                               30.72
                 EAD = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+                      30.72
     &                      SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS*DDIR          30.72
                 ETOT = ETOT + EAD
               ENDDO
               IF (MSC .GT. 3) THEN                                       10.20
!                contribution of tail to total energy density
                 EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                       30.72
                 ETOT = ETOT + DDIR * EHFR * SPCSIG(MSC) * EFTAIL         30.72
               ENDIF
             ENDDO
          ELSE                                                            40.87
!            integration over [fmin,fmax]                                 40.87
             FMIN = PI2*OUTPAR(21)                                        40.87
             FMAX = PI2*OUTPAR(36)                                        40.87
             ECS  = 1.                                                    40.87
             ETOT = SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
          ENDIF                                                           40.87
          IF (ETOT .GE. 0.) THEN                                          30.00
            VOQ(IP,VOQR(IVTYPE)) = 4.*SQRT(ETOT)
          ELSE
            VOQ(IP,VOQR(IVTYPE)) = 0.                                     40.86
          ENDIF
          IF (ITEST.GE.100) THEN                                          40.00
            WRITE(PRINTF, 222) IP, OVSNAM(IVTYPE), VOQ(IP,VOQR(IVTYPE))   40.00
          ENDIF
        ENDIF
!
!       swell wave height
!
        IVTYPE = 44
        IF (OQPROC(IVTYPE)) THEN
          ETOT = 0.
          FSWELL = PI2 * OUTPAR(5)
!         trapezoidal rule is applied
          DO IS = 2, MSC
            DS = SPCSIG(IS)-SPCSIG(IS-1)
            IF (SPCSIG(IS).LE.FSWELL) THEN
              CIA = 0.5 * SPCSIG(IS-1) * DS * DDIR
              CIB = 0.5 * SPCSIG(IS  ) * DS * DDIR
            ELSE
              DSSW = FSWELL-SPCSIG(IS-1)
              CIB = 0.5 * FSWELL * DDIR * DSSW**2 / DS
              CIA = 0.5 * (SPCSIG(IS-1)+FSWELL) * DSSW * DDIR - CIB       40.87
            ENDIF
            DO ID = 1, MDC
              EAD = CIA * ACLOC(ID,IS-1) + CIB * ACLOC(ID,IS)
              ETOT = ETOT + EAD
            ENDDO
            IF (SPCSIG(IS).GT.FSWELL) EXIT
          ENDDO
          IF (ETOT .GE. 0.) THEN                                          30.00
            VOQ(IP,VOQR(IVTYPE)) = 4.*SQRT(ETOT)
          ELSE
            VOQ(IP,VOQR(IVTYPE)) = 0.                                     40.86
          ENDIF
          IF (ITEST.GE.100) THEN                                          40.00
            WRITE(PRINTF, 222) IP, OVSNAM(IVTYPE), VOQ(IP,VOQR(IVTYPE))   40.00
 222        FORMAT(' SWOEXA: POINT ', I5, 2X, A, 1X, E12.4)
          ENDIF
        ENDIF
!
!       average relative period                              modified 10.09
!
        IVTYPE = 28
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(14).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS) ** 2                                40.00
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)                   10.30
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)      30.72
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN                                        10.20
                 PPTAIL = PWTAIL(1) - 1.                                  40.00
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 PPTAIL = PWTAIL(1) - 2.                                  40.00
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)                           10.30
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR                              30.72
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(29)                                       40.87
              FMAX = PI2*OUTPAR(44)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT=SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EPTOT=SwanIntgratSpc(1. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
              TPER = 2.*PI * APTOT / EPTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
        IVTYPE = 11
        IF (ICUR.EQ.0.AND.OQPROC(IVTYPE)) THEN
           IF (OUTPAR(7).EQ.0.) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS) ** 2                                40.00
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)                   10.30
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)      30.72
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN                                        10.20
                 PPTAIL = PWTAIL(1) - 1.                                  40.00
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 PPTAIL = PWTAIL(1) - 2.                                  40.00
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)                           10.30
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR                              30.72
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(22)                                       40.87
              FMAX = PI2*OUTPAR(37)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT=SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EPTOT=SwanIntgratSpc(1. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
             TPER = 2.*PI * APTOT / EPTOT
             VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
             VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       average relative period                              modified 10.09
!
        IVTYPE = 43                                                       40.00
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(18).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS) ** (OUTPAR(2)+1.)                   40.00
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)                   10.30
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)      30.72
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN                                        10.20
                 PPTAIL = PWTAIL(1) - OUTPAR(2)                           40.00
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 PPTAIL = PWTAIL(1) - OUTPAR(2) - 1.                      40.00
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)                           10.30
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR                              30.72
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(33)                                       40.87
              FMAX = PI2*OUTPAR(48)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT = SwanIntgratSpc(OUTPAR(2)-1., FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 1   )           40.87
              EPTOT = SwanIntgratSpc(OUTPAR(2)   , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 1   )           40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
              TPER = 2.*PI * APTOT / EPTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
        IVTYPE = 42                                                       40.00
        IF (ICUR.EQ.0.AND.OQPROC(IVTYPE)) THEN
           IF (OUTPAR(17).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS) ** (OUTPAR(2)+1.)                   40.00
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)                   10.30
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)      30.72
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN                                        10.20
                 PPTAIL = PWTAIL(1) - OUTPAR(2)                           40.00
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 PPTAIL = PWTAIL(1) - OUTPAR(2) - 1.                      40.00
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.61
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)                           10.30
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR                              30.72
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(32)                                       40.87
              FMAX = PI2*OUTPAR(47)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT = SwanIntgratSpc(OUTPAR(2)-1., FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 1   )           40.87
              EPTOT = SwanIntgratSpc(OUTPAR(2)   , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 1   )           40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
              TPER = 2.*PI * APTOT / EPTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       peak period
!
        IVTYPE = 12
        IF (OQPROC(IVTYPE)) THEN
           EMAX = 0.
           ISIGM = -1
           DO IS = 1, MSC
              ETD = 0.
              DO ID = 1, MDC
                ETD = ETD + SPCSIG(IS)*ACLOC(ID,IS)*DDIR                  30.72
              ENDDO
              IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
              ENDIF
           ENDDO
           IF (ISIGM.GT.0) THEN
             VOQ(IP,VOQR(IVTYPE)) = 2.*PI/SPCSIG(ISIGM)                   30.72
           ELSE
             VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       peak period based on parabolic fitting                            40.51
!
        IVTYPE = 53                                                       40.51
        IF (OQPROC(IVTYPE)) THEN
           EMAX = 0.
           ETD  = 0.
           ISIGM = -1
           DO IS = 1, MSC
              ED  = ETD
              ETD = 0.
              DO ID = 1, MDC
                ETD = ETD + SPCSIG(IS)*ACLOC(ID,IS)*DDIR
              END DO
              IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
                EMAXD = ED
                EMAXU = 0.
                IF (IS.LT.MSC) THEN
                   DO ID = 1, MDC
                     EMAXU = EMAXU + SPCSIG(IS+1)*ACLOC(ID,IS+1)*DDIR
                   END DO
                ELSE
                   EMAXU = EMAX
                END IF
              END IF
           END DO
           IF (ISIGM.GT.1 .AND. ISIGM.LT.MSC) THEN                        41.20
             SIG1 = SPCSIG(ISIGM-1)
             SIG2 = SPCSIG(ISIGM+1)
             SIG3 = SPCSIG(ISIGM  )
             E1   = EMAXD
             E2   = EMAXU
             E3   = EMAX
             P    = SIG1+SIG2
             Q    = (E1-E2)/(SIG1-SIG2)
             R    = SIG1+SIG3
             T    = (E1-E3)/(SIG1-SIG3)
             A    = (T-Q)/(R-P)
             IF (A.LT.0) THEN
                SIGP = (-Q+P*A)/(2.*A)
             ELSE
                SIGP = SIG3
             END IF
             VOQ(IP,VOQR(IVTYPE)) = 2.*PI/SIGP
           ELSE IF (ISIGM.EQ.1) THEN
             VOQ(IP,VOQR(IVTYPE)) = 2.*PI/SPCSIG(1)
           ELSE
             VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           END IF
        END IF
!
!       peak direction
!
        IVTYPE = 14
        IF (OQPROC(IVTYPE)) THEN
           EMAX = 0.
           IDIRM = -1
           DO ID = 1, MDC
              ETF = 0.
              DO IS = 2, MSC
                DS = SPCSIG(IS)-SPCSIG(IS-1)                              30.72
                E1 = SPCSIG(IS-1)*ACLOC(ID,IS-1)                          30.72
                E2 = SPCSIG(IS)*ACLOC(ID,IS)                              30.72
                ETF = ETF + DS * (E1+E2)
              ENDDO
              IF (ETF.GT.EMAX) THEN
                EMAX  = ETF
                IDIRM = ID
              ENDIF
           ENDDO
           IF (IDIRM.GT.0) THEN
!
!            *** Convert (if necessary) from nautical degrees ***         32.01
!            *** to cartesian degrees                         ***         32.01
!
             IF (BNAUT) THEN                                              30.70
               VOQ(IP,VOQR(IVTYPE)) = ANGDEG( SPCDIR(IDIRM,1) )           32.01
             ELSE
               VOQ(IP,VOQR(IVTYPE)) = ANGDEG( (ALCQ + SPCDIR(IDIRM,1)) )  32.01
             ENDIF
             VOQ(IP,VOQR(IVTYPE)) = DEGCNV( VOQ(IP,VOQR(IVTYPE)) )        32.01
!
           ELSE
             VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       mean direction
!
        IVTYPE = 13
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(8).EQ.0.) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              ETOT = 0.
              EEX  = 0.
              EEY  = 0.
              DO ID=1, MDC
                 EAD = 0.
                 DO IS=2,MSC
                   DS=SPCSIG(IS)-SPCSIG(IS-1)                             30.72
                   EDI = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+                    30.72
     &                        SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS             30.72
                   EAD = EAD + EDI
                 ENDDO
                 IF (MSC .GT. 3) THEN                                     10.20
!                  contribution of tail to total energy density
                   EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                     30.72
                   EAD = EAD + EHFR * SPCSIG(MSC) * EFTAIL                30.72
                 ENDIF
                 EAD = EAD * DDIR
                 ETOT = ETOT + EAD
                 EEX  = EEX + EAD * SPCDIR(ID,2)
                 EEY  = EEY + EAD * SPCDIR(ID,3)
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(23)                                       40.87
              FMAX = PI2*OUTPAR(38)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT= SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EEX = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             WK, SPCDIR(1,2), 0., 0., ACLOC, 1)     40.87
              EEY = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             WK, SPCDIR(1,3), 0., 0., ACLOC, 1)     40.87
           ENDIF                                                          40.87
           IF (ETOT.GT.0.) THEN
              IF (BNAUT) THEN                                             30.70
                 DIRDEG = ATAN2(EEY,EEX) * 180./PI                        10.15
              ELSE
                 DIRDEG = (ALCQ + ATAN2(EEY,EEX)) * 180./PI               10.15
              ENDIF
              IF (DIRDEG.LT.0.) DIRDEG = DIRDEG + 360.                    10.15
!
!             *** Convert (if necessary) from nautical degrees ***        32.01
!             *** to cartesian degrees                         ***        32.01
!
              VOQ(IP,VOQR(IVTYPE)) = DEGCNV( DIRDEG )                     32.01
!
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       directional spread
!
        IVTYPE = 16
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(10).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              ETOT = 0.
              EEX  = 0.
              EEY  = 0.
              DO ID=1, MDC
                 EAD = 0.
                 DO IS=2,MSC
                   DS=SPCSIG(IS)-SPCSIG(IS-1)                             30.72
                   EDI = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+                    30.72
     &                        SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS             30.72
                   EAD = EAD + EDI
                 ENDDO
                 IF (MSC .GT. 3) THEN                                     10.20
!                  contribution of tail to total energy density
                   EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                     30.72
                   EAD = EAD + EHFR * SPCSIG(MSC) * EFTAIL                30.72
                 ENDIF
                 EAD = EAD * DDIR
                 ETOT = ETOT + EAD
                 EEX  = EEX + EAD * SPCDIR(ID,2)
                 EEY  = EEY + EAD * SPCDIR(ID,3)
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(25)                                       40.87
              FMAX = PI2*OUTPAR(40)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT= SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EEX = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             WK, SPCDIR(1,2), 0., 0., ACLOC, 1)     40.87
              EEY = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             WK, SPCDIR(1,3), 0., 0., ACLOC, 1)     40.87
           ENDIF                                                          40.87
           IF (ETOT.GT.0.) THEN
              FF = MIN (1., SQRT(EEX*EEX+EEY*EEY)/ETOT)
              VOQ(IP,VOQR(IVTYPE)) = SQRT(2.-2.*FF) *180./PI
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       average relative period (Tm-10)                                   40.41
!
        IVTYPE = 48
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(20).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS)
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN
                 PPTAIL = PWTAIL(1)
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
                 PPTAIL = PWTAIL(1) - 1.
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(35)                                       40.87
              FMAX = PI2*OUTPAR(50)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT=SwanIntgratSpc(-1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EPTOT=SwanIntgratSpc( 0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
              TPER = 2.*PI * APTOT / EPTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
        IVTYPE = 47
        IF (ICUR.EQ.0.AND.OQPROC(IVTYPE)) THEN
           IF (OUTPAR(19).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              APTOT = 0.
              EPTOT = 0.
              DO ID=1, MDC
                 DO IS=1,MSC
                   SIG2P = SPCSIG(IS)
                   APTOT = APTOT + SIG2P * ACLOC(ID,IS)
                   EPTOT = EPTOT + SPCSIG(IS) * SIG2P * ACLOC(ID,IS)
                 ENDDO
              ENDDO
              APTOT = APTOT * FRINTF
              EPTOT = EPTOT * FRINTF
              IF (MSC .GT. 3) THEN
                 PPTAIL = PWTAIL(1)
                 APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
                 PPTAIL = PWTAIL(1) - 1.
                 EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
                 DO ID = 1, MDC
!                  contribution of tail to total energy density
                   AHFR = SIG2P * ACLOC(ID,MSC)
                   APTOT = APTOT + APTAIL * AHFR
                   EHFR = SPCSIG(MSC) * AHFR
                   EPTOT = EPTOT + EPTAIL * EHFR
                 ENDDO
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(34)                                       40.87
              FMAX = PI2*OUTPAR(49)                                       40.87
              ECS  = 1.                                                   40.87
              APTOT=SwanIntgratSpc(-1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
              EPTOT=SwanIntgratSpc( 0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
     &                             1  )                                   40.87
           ENDIF                                                          40.87
           IF (EPTOT.GT.0.) THEN
              TPER = 2.*PI * APTOT / EPTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       peakedness (Qp)                                                   40.64
!
        IVTYPE = 58
        IF (OQPROC(IVTYPE)) THEN
           ETOT  = 0.
           ESTOT = 0.
           DO ID=1, MDC
             DO IS = 1, MSC
               SIG   = SPCSIG(IS)
               EADD  = SIG**2 * ACLOC(ID,IS) * FRINTF * DDIR
               ETOT  = ETOT  + EADD
             ENDDO
           ENDDO
           DO IS = 1, MSC
             SIG = SPCSIG(IS)
             EADD = 0.
             DO ID=1, MDC
               EADD = EADD + SIG * ACLOC(ID,IS) * DDIR
             ENDDO
             ESTOT = ESTOT + SIG**2 * EADD**2 * FRINTF
           ENDDO
           IF (ETOT.GT.0.) THEN
              VOQ(IP,VOQR(IVTYPE)) = 2.*ESTOT/(ETOT * ETOT)
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
        IF (BKC.EQ.1) GOTO 800
!
!       compute k and Cg
!
        DEPLOC = VOQ(IP,VOQR(4))
        CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)                30.81 30.72
        IF (ITEST.GE.100 .OR. IOUTES .GE. 20) THEN
          WRITE (PRTEST, *) ' Depth: ', DEPLOC
          DO 870 ISC = 1, MIN(MSC,20)
             WRITE (PRTEST, 860) ISC, SPCSIG(ISC), WK(ISC), CG(ISC)       30.72
 860         FORMAT (' i, SPCSIG, sigma, k, cg ', I2, 3(1X, E12.4))
 870      CONTINUE
        ENDIF
!
!       transport direction
!
        IVTYPE = 15
        IF (OQPROC(IVTYPE)) THEN
           IF (ICUR.EQ.0) THEN
             UXLOC = 0.
             UYLOC = 0.
           ELSE
             UXLOC = VOQ(IP,VOQR(5))
             UYLOC = VOQ(IP,VOQR(5)+1)
           ENDIF
           IF (OUTPAR(9).EQ.0) THEN                                       40.87
!             integration over [0,inf]                                    40.87
              CEX = 0.
              CEY = 0.
              ETOT = 0.
              DO ISIGM = 1, MSC
                IF (ISIGM.EQ.1) THEN
                  DSIG = 0.5 * (SPCSIG(2) - SPCSIG(1))                    30.72
                ELSE IF (ISIGM.EQ.MSC) THEN
                  DSIG = 0.5 * (SPCSIG(MSC) - SPCSIG(MSC-1))              30.72
                ELSE
                  DSIG = 0.5 * (SPCSIG(ISIGM+1) - SPCSIG(ISIGM-1))        30.72
                ENDIF
                SIG2 = SPCSIG(ISIGM)                                      30.72
                CS   = CG(ISIGM)*SIG2
                DO ID=1,MDC
                   CGE = DSIG * CS * ACLOC(ID,ISIGM)
                   CEX = CEX + CGE * SPCDIR(ID,2)
                   CEY = CEY + CGE * SPCDIR(ID,3)
                   IF (ICUR.EQ.1) THEN
                     ETOT = ETOT + DSIG * SIG2 * ACLOC(ID,ISIGM)
                   ENDIF
                ENDDO
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(24)                                       40.87
              FMAX = PI2*OUTPAR(39)                                       40.87
              ECS  = 1.                                                   40.87
              IF (ICUR.EQ.1)                                              40.87
     &           ETOT=SwanIntgratSpc(0.,FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                               WK, ECS, 0., 0., ACLOC, 1)           40.87
              CEX = SwanIntgratSpc(1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             CG, SPCDIR(1,2), 0., 0., ACLOC, 4)     40.87
              CEY = SwanIntgratSpc(1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             CG, SPCDIR(1,3), 0., 0., ACLOC, 4)     40.87
           ENDIF                                                          40.87
!
           IF (ICUR.EQ.1) THEN
              CEX = CEX + ETOT * UXLOC
              CEY = CEY + ETOT * UYLOC
           ENDIF
!
           IF (OQPROC(IVTYPE)) THEN
              IF (CEX.EQ.0. .AND. CEY.EQ.0.) THEN
                VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
              ELSE
                IF (BNAUT) THEN                                           30.70
                  DIRDEG = ATAN2(CEY,CEX) * 180./PI                       10.15
                ELSE
                  DIRDEG = (ALCQ + ATAN2(CEY,CEX)) * 180./PI              10.15
                ENDIF
                IF (DIRDEG.LT.0.) DIRDEG = DIRDEG + 360.                  10.15
!
!               *** Convert (if necessary) from nautical degrees ***      32.01
!               *** to cartesian degrees                         ***      32.01
!
                VOQ(IP,VOQR(IVTYPE)) = DEGCNV( DIRDEG )                   32.01
!
              ENDIF
           ENDIF
        ENDIF
!
!       transport vector
!
        IVTYPE = 19
        IF (OQPROC(IVTYPE)) THEN
           IF (ICUR.EQ.0) THEN
             UXLOC = 0.
             UYLOC = 0.
           ELSE
             UXLOC = VOQ(IP,VOQR(5))
             UYLOC = VOQ(IP,VOQR(5)+1)
           ENDIF
           IF (OUTPAR(13).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              CEX = 0.
              CEY = 0.
              ETOT = 0.
              DO ISIGM = 1, MSC
                IF (ISIGM.EQ.1) THEN
                  DSIG = 0.5 * (SPCSIG(2) - SPCSIG(1))                    30.72
                ELSE IF (ISIGM.EQ.MSC) THEN
                  DSIG = 0.5 * (SPCSIG(MSC) - SPCSIG(MSC-1))              30.72
                ELSE
                  DSIG = 0.5 * (SPCSIG(ISIGM+1) - SPCSIG(ISIGM-1))        30.72
                ENDIF
                SIG2 = SPCSIG(ISIGM)                                      30.72
                CS   = CG(ISIGM)*SIG2
                DO ID=1,MDC
                   CGE = DSIG * CS * ACLOC(ID,ISIGM)
                   CEX = CEX + CGE * SPCDIR(ID,2)
                   CEY = CEY + CGE * SPCDIR(ID,3)
                   IF (ICUR.EQ.1) THEN
                     ETOT = ETOT + DSIG * SIG2 * ACLOC(ID,ISIGM)
                   ENDIF
                ENDDO
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(28)                                       40.87
              FMAX = PI2*OUTPAR(43)                                       40.87
              ECS  = 1.                                                   40.87
              IF (ICUR.EQ.1)                                              40.87
     &           ETOT=SwanIntgratSpc(0.,FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                               WK, ECS, 0., 0., ACLOC, 1)           40.87
              CEX = SwanIntgratSpc(1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             CG, SPCDIR(1,2), 0., 0., ACLOC, 4)     40.87
              CEY = SwanIntgratSpc(1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),   40.87
     &                             CG, SPCDIR(1,3), 0., 0., ACLOC, 4)     40.87
              CEX = CEX/DDIR                                              40.87
              CEY = CEY/DDIR                                              40.87
           ENDIF                                                          40.87
!
           IF (ICUR.EQ.1) THEN
              CEX = CEX + ETOT * UXLOC
              CEY = CEY + ETOT * UYLOC
           ENDIF
!
           SX = CEX * DDIR
           SY = CEY * DDIR
           IF (INRHOG.EQ.1) THEN
              SX = SX * RHO * GRAV
              SY = SY * RHO * GRAV
           ENDIF
           VOQ(IP,VOQR(IVTYPE))   = COSCQ*SX - SINCQ*SY
           VOQ(IP,VOQR(IVTYPE)+1) = SINCQ*SX + COSCQ*SY
        ENDIF
!
!       average wave length
!
        IVTYPE = 17
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(11).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              ETOT  = 0.
              EKTOT = 0.
!             new integration method involving FRINTF                     20.59
              DO IS=1, MSC
                 SIG2 = (SPCSIG(IS))**2                                   30.72
                 SKK  = SIG2 * (WK(IS))**OUTPAR(3)                        40.00
                 DO ID=1,MDC
                   ETOT  = ETOT + SIG2 * ACLOC(ID,IS)                     20.59
                   EKTOT = EKTOT + SKK * ACLOC(ID,IS)                     20.59
                 ENDDO
              ENDDO
              ETOT  = FRINTF * ETOT
              EKTOT = FRINTF * EKTOT
              IF (MSC .GT. 3) THEN                                        10.20
!                contribution of tail to total energy density
                 PPTAIL = PWTAIL(1) - 1.                                  20.59
                 CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.59
                 PPTAIL = PWTAIL(1) - 1. - 2.*OUTPAR(3)                   40.00
                 IF (PPTAIL.LE.0.) THEN
                   CALL MSGERR (2,'error tail computation')
                   GOTO 480
                 ENDIF
                 CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.59
                 DO ID=1,MDC
                   ETOT   = ETOT + CETAIL * SIG2 * ACLOC(ID,MSC)          20.59
                   EKTOT  = EKTOT + CKTAIL * SKK * ACLOC(ID,MSC)          20.59
                 ENDDO
 480             CONTINUE
              ENDIF
              IF (EKTOT.GT.0.) THEN
                 WLMEAN = PI2 * (ETOT / EKTOT) ** (1./OUTPAR(3))          40.00
                 VOQ(IP,VOQR(IVTYPE)) = WLMEAN
              ELSE
                 VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(26)                                       40.87
              FMAX = PI2*OUTPAR(41)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT  = SwanIntgratSpc(OUTPAR(3)-1., FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 3   )           40.87
              EKTOT = SwanIntgratSpc(OUTPAR(3)   , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 3   )           40.87
              IF (EKTOT.GT.0.) THEN                                       40.87
                 WLMEAN = PI2 * ETOT / EKTOT                              40.87
                 VOQ(IP,VOQR(IVTYPE)) = WLMEAN                            40.87
              ELSE                                                        40.87
                 VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)                    40.87
              ENDIF                                                       40.87
           ENDIF                                                          40.87
        ENDIF
!
!       steepness
!
        IVTYPE = 18
        IF (OQPROC(IVTYPE)) THEN
           IF (OUTPAR(12).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              ETOT  = 0.
              EKTOT = 0.
!             new integration method involving FRINTF                     20.59
              DO IS=1, MSC
                 SIG2 = (SPCSIG(IS))**2                                   30.72
                 SKK  = SIG2 * (WK(IS))**OUTPAR(3)                        40.00
                 DO ID=1,MDC
                   ETOT  = ETOT + SIG2 * ACLOC(ID,IS)                     20.59
                   EKTOT = EKTOT + SKK * ACLOC(ID,IS)                     20.59
                 ENDDO
              ENDDO
              ETOT  = FRINTF * ETOT
              EKTOT = FRINTF * EKTOT
              IF (MSC .GT. 3) THEN                                        10.20
!                contribution of tail to total energy density
                 PPTAIL = PWTAIL(1) - 1.                                  20.59
                 CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.59
                 PPTAIL = PWTAIL(1) - 1. - 2.*OUTPAR(3)                   40.00
                 IF (PPTAIL.LE.0.) THEN
                   CALL MSGERR (2,'error tail computation')
                   GOTO 481
                 ENDIF
                 CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))     20.59
                 DO ID=1,MDC
                   ETOT   = ETOT + CETAIL * SIG2 * ACLOC(ID,MSC)          20.59
                   EKTOT  = EKTOT + CKTAIL * SKK * ACLOC(ID,MSC)          20.59
                 ENDDO
 481             CONTINUE
              ENDIF
              IF (EKTOT.GT.0.) THEN
                 WLMEAN = PI2 * (ETOT / EKTOT) ** (1./OUTPAR(3))          40.00
                 VOQ(IP,VOQR(IVTYPE)) = 4.* SQRT(ETOT*DDIR) / WLMEAN
              ELSE
                 VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
              ENDIF
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(27)                                       40.87
              FMAX = PI2*OUTPAR(42)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT  = SwanIntgratSpc(0.          , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 1   )           40.87
              EPTOT = SwanIntgratSpc(OUTPAR(3)-1., FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 3   )           40.87
              EKTOT = SwanIntgratSpc(OUTPAR(3)   , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , 0.    ,   40.87
     &                               0.          , ACLOC, 3   )           40.87
              IF (EKTOT.GT.0.) THEN                                       40.87
                 WLMEAN = PI2 * EPTOT / EKTOT                             40.87
                 VOQ(IP,VOQR(IVTYPE)) = 4.* SQRT(ETOT) / WLMEAN           40.87
              ELSE                                                        40.87
                 VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)                    40.87
              ENDIF                                                       40.87
           ENDIF                                                          40.87
        ENDIF
!
!       average absolute period Tm01                                      40.00
!
        IVTYPE = 11
        IF (ICUR.GT.0 .AND. OQPROC(IVTYPE)) THEN
           UXLOC = VOQ(IP,VOQR(5))
           UYLOC = VOQ(IP,VOQR(5)+1)
           IF (OUTPAR(7).EQ.0) THEN                                       40.87
!             integration over [0,inf]                                    40.87
              ETOT = 0.
              EFTOT = 0.
              PPTAIL = PWTAIL(1) - 1.                                     40.00
              ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              PPTAIL = PWTAIL(1) - 2.                                     40.00
              EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              DO ID=1, MDC
                 THETA = SPCDIR(ID,1) + ALCQ                              20.43
                 UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 DO IS = 1, MSC
                   OMEG = SPCSIG(IS) + WK(IS) * UXD                       30.72
                   EADD = FRINTF * SPCSIG(IS)**2 * ACLOC(ID,IS)           40.00
                   ETOT = ETOT + EADD
                   EFTOT = EFTOT + EADD * OMEG                            20.66
                 ENDDO
                 IF (MSC .GT. 3) THEN                                     10.20
!                  contribution of tail to total energy density
                   EADD = SPCSIG(MSC)**2 * ACLOC(ID,MSC)                  40.00
                   ETOT = ETOT + ETAIL * EADD
                   EFTOT = EFTOT + EFTAIL * OMEG * EADD
                 ENDIF
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(22)                                       40.87
              FMAX = PI2*OUTPAR(37)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT =SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , UXLOC, UYLOC, ACLOC      ,  40.87
     &                             2  )                                   40.87
              EFTOT=SwanIntgratSpc(1. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , UXLOC, UYLOC, ACLOC      ,  40.87
     &                             2  )                                   40.87
           ENDIF                                                          40.87
           IF (EFTOT.GT.0.) THEN
              TPER = 2.*PI * ETOT / EFTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       average absolute period (case with current)                       10.09
!
        IVTYPE = 42                                                       40.00
        IF (ICUR.GT.0 .AND. OQPROC(IVTYPE)) THEN
           UXLOC = VOQ(IP,VOQR(5))
           UYLOC = VOQ(IP,VOQR(5)+1)
           IF (OUTPAR(17).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              ETOT = 0.
              EFTOT = 0.
              PPTAIL = PWTAIL(1) - OUTPAR(2)                              40.00
              ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              PPTAIL = PWTAIL(1) - OUTPAR(2) - 1.                         40.00
              EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              DO ID=1, MDC
                 THETA = SPCDIR(ID,1) + ALCQ                              20.43
                 UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 DO IS = 1, MSC
                   OMEG = SPCSIG(IS) + WK(IS) * UXD                       10.30
                   OMEG1P = OMEG ** (OUTPAR(2)-1.)                        10.30
                   EADD = OMEG1P * FRINTF * SPCSIG(IS)**2 * ACLOC(ID,IS)  20.66
                   ETOT = ETOT + EADD
                   EFTOT = EFTOT + EADD * OMEG                            20.66
                 ENDDO
                 IF (MSC .GT. 3) THEN                                     10.20
!                  contribution of tail to total energy density
                   EADD = OMEG1P * SPCSIG(MSC)**2 * ACLOC(ID,MSC)
                   ETOT = ETOT + ETAIL * EADD
                   EFTOT = EFTOT + EFTAIL * OMEG * EADD
                 ENDIF
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(32)                                       40.87
              FMAX = PI2*OUTPAR(47)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT  = SwanIntgratSpc(OUTPAR(2)-1., FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , UXLOC ,   40.87
     &                               UYLOC       , ACLOC, 2   )           40.87
              EFTOT = SwanIntgratSpc(OUTPAR(2)   , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , UXLOC ,   40.87
     &                               UYLOC       , ACLOC, 2   )           40.87
           ENDIF                                                          40.87
           IF (EFTOT.GT.0.) THEN
              TPER = 2.*PI * ETOT / EFTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       zero-crossing period Tm02                                         20.61
!
        IVTYPE = 32                                                       20.61
        IF (OQPROC(IVTYPE)) THEN
           IF (ICUR.GT.0) THEN
             UXLOC = VOQ(IP,VOQR(5))
             UYLOC = VOQ(IP,VOQR(5)+1)
           ENDIF
           IF (OUTPAR(15).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              ETOT  = 0.
              EFTOT = 0.
              PPTAIL = PWTAIL(1) - 1.                                     20.61
              ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              PPTAIL = PWTAIL(1) - 3.                                     20.61
              EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        20.61
              DO ID=1, MDC
                 IF (ICUR.GT.0) THEN
                   THETA = SPCDIR(ID,1) + ALCQ
                   UXD   = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 ENDIF
                 DO IS=1,MSC
                   EADD  = SPCSIG(IS)**2 * ACLOC(ID,IS) * FRINTF          30.72
                   IF (ICUR.GT.0) THEN
                     OMEG  = SPCSIG(IS) + WK(IS) * UXD                    30.72
                     OMEG2 = OMEG**2
                   ELSE
                     OMEG2 = SPCSIG(IS)**2                                30.72
                   ENDIF
                   ETOT  = ETOT + EADD                                    20.61
                   EFTOT = EFTOT + EADD * OMEG2                           20.61
                 ENDDO
                 IF (MSC .GT. 3) THEN
!                  contribution of tail to total energy density
                   EADD  = SPCSIG(MSC)**2 * ACLOC(ID,MSC)                 30.72
                   ETOT  = ETOT  + ETAIL * EADD                           20.61
                   EFTOT = EFTOT + EFTAIL * OMEG2 * EADD                  20.61
                 ENDIF
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(30)                                       40.87
              FMAX = PI2*OUTPAR(45)                                       40.87
              ECS  = 1.                                                   40.87
              IF (ICUR.GT.0) THEN                                         40.87
                 ITP = 2                                                  40.87
              ELSE                                                        40.87
                 ITP = 1                                                  40.87
              ENDIF                                                       40.87
              ETOT  = SwanIntgratSpc(0.          , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , UXLOC ,   40.87
     &                               UYLOC       , ACLOC, ITP )           40.87
              EFTOT = SwanIntgratSpc(2.          , FMIN , FMAX, SPCSIG,   40.87
     &                               SPCDIR(1,1) , WK   , ECS , UXLOC ,   40.87
     &                               UYLOC       , ACLOC, ITP )           40.87
           ENDIF                                                          40.87
           IF (EFTOT.GT.0.) THEN
              VOQ(IP,VOQR(IVTYPE)) = 2.*PI * SQRT(ETOT/EFTOT)             20.61
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       frequency spectral width (kappa)                                  20.61
!
        IVTYPE = 33                                                       20.61
        IF (OQPROC(IVTYPE)) THEN
           TM02 = VOQ(IP,VOQR(32))
           IF (ICUR.GT.0) THEN
             UXLOC = VOQ(IP,VOQR(5))
             UYLOC = VOQ(IP,VOQR(5)+1)
           ENDIF
           ETOT  = 0.
           ECTOT = 0.
           ESTOT = 0.
           IF (OUTPAR(16).EQ.0) THEN                                      40.87
!             integration over [0,inf]                                    40.87
              PPTAIL = PWTAIL(1) - 1.
              ECTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
              DO  ID=1, MDC
                IF (ICUR.GT.0) THEN
                  THETA = SPCDIR(ID,1) + ALCQ
                  UXD   = UXLOC*COS(THETA) + UYLOC*SIN(THETA)             20.66
                ENDIF
                DO  IS = 1, MSC
                  SIG = SPCSIG(IS)                                        30.72
                  IF (ICUR.GT.0) THEN
                    OMEG = SIG + WK(IS) * UXD                             20.66
                  ELSE
                    OMEG = SIG
                  ENDIF
                  FND = OMEG * TM02
                  COSFND = COS(FND)
                  SINFND = SIN(FND)
                  EADD   = SIG**2 * ACLOC(ID,IS) * FRINTF                 30.70
                  ETOT  = ETOT  + EADD
                  ECTOT = ECTOT + COSFND * EADD                           20.66
                  ESTOT = ESTOT + SINFND * EADD                           20.66
                ENDDO
                IF (MSC .GT. 3) THEN
!                 contribution of tail to total energy density
                  EADD  = ECTAIL * SIG**2 * ACLOC(ID,MSC)                 30.70
                  ETOT  = ETOT  + EADD
                  ECTOT = ECTOT + COSFND * EADD                           20.66
                  ESTOT = ESTOT + SINFND * EADD                           20.66
                ENDIF
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(31)                                       40.87
              FMAX = PI2*OUTPAR(46)                                       40.87
              DO ID = 1, MDC
                IF (ICUR.GT.0) THEN
                  THETA = SPCDIR(ID,1) + ALCQ
                  UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                ENDIF
                DO IS = 2, MSC
                  SIG1 = SPCSIG(IS-1)
                  SIG2 = SPCSIG(IS  )
                  IF (ICUR.GT.0) THEN
                    OMEG1 = SIG1 + WK(IS-1) * UXD
                    OMEG2 = SIG2 + WK(IS  ) * UXD
                  ELSE
                    OMEG1 = SIG1
                    OMEG2 = SIG2
                  ENDIF
                  COSFN1 = COS(OMEG1 * TM02)
                  SINFN1 = SIN(OMEG1 * TM02)
                  COSFN2 = COS(OMEG2 * TM02)
                  SINFN2 = SIN(OMEG2 * TM02)
                  DS = SIG2 - SIG1
                  IF ( SIG1.GE.FMIN .AND. SIG2.LE.FMAX ) THEN
                    CIA  = 0.5*SIG1 * DS
                    CIB  = 0.5*SIG2 * DS
                    CIAC = 0.5*COSFN1 * SIG1 * DS
                    CIBC = 0.5*COSFN2 * SIG2 * DS
                    CIAS = 0.5*SINFN1 * SIG1 * DS
                    CIBS = 0.5*SINFN2 * SIG2 * DS
                  ELSEIF ( SIG2.GT.FMAX ) THEN
                    DSSW = FMAX - SIG1
                    IF (ICUR.GT.0) THEN
                      OMEG2=FMAX+(WK(IS)*DSSW+WK(IS-1)*(DS-DSSW))*UXD/DS
                    ELSE
                      OMEG2=FMAX
                    ENDIF
                    COSFN2 = COS(OMEG2 * TM02)
                    SINFN2 = SIN(OMEG2 * TM02)
                    CIBC   = 0.5*COSFN2 * FMAX * DSSW**2 / DS
                    CIAC   = 0.5*(COSFN1*SIG1 + COSFN2*FMAX)*DSSW - CIBC
                    CIBS   = 0.5*SINFN2 * FMAX * DSSW**2 / DS
                    CIAS   = 0.5*(SINFN1*SIG1 + SINFN2*FMAX)*DSSW - CIBS
                    CIB    = 0.5*FMAX * DSSW**2 / DS
                    CIA    = 0.5*(SIG1 + FMAX)*DSSW - CIB
                  ELSEIF ( SIG2.GT.FMIN ) THEN
                    DSSW = SIG2 - FMIN
                    IF (ICUR.GT.0) THEN
                      OMEG1=FMIN+(WK(IS-1)*DSSW+WK(IS)*(DS-DSSW))*UXD/DS
                    ELSE
                      OMEG1=FMIN
                    ENDIF
                    COSFN1 = COS(OMEG1 * TM02)
                    SINFN1 = SIN(OMEG1 * TM02)
                    CIAC   = 0.5*COSFN1 * FMIN * DSSW**2 / DS
                    CIBC   = 0.5*(COSFN2*SIG2 + COSFN1*FMIN)*DSSW - CIAC
                    CIAS   = 0.5*SINFN1 * FMIN * DSSW**2 / DS
                    CIBS   = 0.5*(SINFN2*SIG2 + SINFN1*FMIN)*DSSW - CIAS
                    CIA    = 0.5*FMIN * DSSW**2 / DS
                    CIB    = 0.5*(SIG2 + FMIN)*DSSW - CIA
                  ELSE
                    CIA  = 0.
                    CIB  = 0.
                    CIAC = 0.
                    CIBC = 0.
                    CIAS = 0.
                    CIBS = 0.
                  ENDIF
                  ETOT = ETOT  + CIA *ACLOC(ID,IS-1) + CIB *ACLOC(ID,IS)
                  ECTOT= ECTOT + CIAC*ACLOC(ID,IS-1) + CIBC*ACLOC(ID,IS)
                  ESTOT= ESTOT + CIAS*ACLOC(ID,IS-1) + CIBS*ACLOC(ID,IS)
                  IF ( SIG2.GT.FMAX ) EXIT
                ENDDO
              ENDDO
!             --- add tail contribution, if appropriate
              IF ( FMAX.GT.SPCSIG(MSC) ) THEN
                 IF ( MSC.GT.3 ) THEN
                    ECTAIL = 1. / (PWTAIL(1) - 1.)
                    IF ( FMAX.GT.100. ) THEN
                       CTAIL = 0.
                    ELSE
                       CTAIL = FMAX * (SPCSIG(MSC)/FMAX)**PWTAIL(1)
                    ENDIF
                    DO ID = 1, MDC
                       IF (ICUR.GT.0) THEN
                          THETA = SPCDIR(ID,1) + ALCQ
                          UXD   = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                          OMEG2 = SPCSIG(MSC) + WK(MSC) * UXD
                       ELSE
                          OMEG2 = SPCSIG(MSC)
                       ENDIF
                       COSFN2 = COS(OMEG2 * TM02)
                       SINFN2 = SIN(OMEG2 * TM02)
                       EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)
                       ETOT  = ETOT  + EHFR*(SPCSIG(MSC)-CTAIL) * ECTAIL
                       ECTOT = ECTOT + EHFR*(COSFN2*SPCSIG(MSC) -
     &                                      COS(FMAX*TM02)*CTAIL)*ECTAIL
                       ESTOT = ESTOT + EHFR*(SINFN2*SPCSIG(MSC) -
     &                                      SIN(FMAX*TM02)*CTAIL)*ECTAIL
                    ENDDO
                 ENDIF
              ENDIF
           ENDIF                                                          40.87
           IF (ETOT.GT.0.) THEN
              VOQ(IP,VOQR(IVTYPE)) =
     &                       SQRT(ECTOT*ECTOT+ESTOT*ESTOT) / ETOT
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       average absolute period Tm-10 (case with current)
!
        IVTYPE = 47
        IF (ICUR.GT.0 .AND. OQPROC(IVTYPE)) THEN
           UXLOC = VOQ(IP,VOQR(5))
           UYLOC = VOQ(IP,VOQR(5)+1)
           IF (OUTPAR(19).EQ.0.) THEN                                     40.87
!             integration over [0,inf]                                    40.87
              ETOT = 0.
              EFTOT = 0.
              PPTAIL = PWTAIL(1)
              ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
              PPTAIL = PWTAIL(1) - 1.
              EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
              DO ID=1, MDC
                 THETA = SPCDIR(ID,1) + ALCQ
                 UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 DO IS = 1, MSC
                   OMEG = SPCSIG(IS) + WK(IS) * UXD
                   OMEG1P = OMEG ** (-1.)
                   EADD = OMEG1P * FRINTF * SPCSIG(IS)**2 * ACLOC(ID,IS)
                   ETOT = ETOT + EADD
                   EFTOT = EFTOT + EADD * OMEG
                 ENDDO
                 IF (MSC .GT. 3) THEN
!                  contribution of tail to total energy density
                   EADD = OMEG1P * SPCSIG(MSC)**2 * ACLOC(ID,MSC)
                   ETOT = ETOT + ETAIL * EADD
                   EFTOT = EFTOT + EFTAIL * OMEG * EADD
                 ENDIF
              ENDDO
           ELSE                                                           40.87
!             integration over [fmin,fmax]                                40.87
              FMIN = PI2*OUTPAR(34)                                       40.87
              FMAX = PI2*OUTPAR(49)                                       40.87
              ECS  = 1.                                                   40.87
              ETOT =SwanIntgratSpc(-1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , UXLOC, UYLOC, ACLOC      ,  40.87
     &                             2  )                                   40.87
              EFTOT=SwanIntgratSpc( 0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
     &                             WK , ECS , UXLOC, UYLOC, ACLOC      ,  40.87
     &                             2  )                                   40.87
           ENDIF                                                          40.87
           IF (EFTOT.GT.0.) THEN
              TPER = 2.*PI * ETOT / EFTOT
              VOQ(IP,VOQR(IVTYPE)) = TPER
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       Benjamin-Feir Index (BFI)                                         40.64
!
        IVTYPE = 59
        IF (OQPROC(IVTYPE)) THEN
           STPNS = VOQ(IP,VOQR(18))
           QP    = VOQ(IP,VOQR(58))
           IF ( STPNS.NE.OVEXCV(18) .AND. QP.NE.OVEXCV(58) ) THEN
              VOQ(IP,VOQR(IVTYPE)) = SQRT(PI2)*STPNS*QP
           ELSE
              VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
!       peak wave length                                                  41.15
!
        IVTYPE = 71
        IF (OQPROC(IVTYPE)) THEN
           EMAX = 0.
           ISIGM = -1
           DO IS = 1, MSC
              ETD = 0.
              DO ID = 1, MDC
                ETD = ETD + WK(IS)*ACLOC(ID,IS)*DDIR
              ENDDO
              IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
              ENDIF
           ENDDO
           IF (ISIGM.GT.0) THEN
             VOQ(IP,VOQR(IVTYPE)) = 2.*PI/WK(ISIGM)
           ELSE
             VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
           ENDIF
        ENDIF
!
        GOTO 800
!
!       points on land: assign exception value
!
 700    DO 730 II = 1, NVOTP
          IVTYPE = IVOTP(II)
          IF (OQPROC(IVTYPE)) THEN
            VOQ(IP,VOQR(IVTYPE)) = OVEXCV(IVTYPE)
            IF (OVSVTY(IVTYPE).EQ.3) THEN
              VOQ(IP,VOQR(IVTYPE)+1) = OVEXCV(IVTYPE)
            ENDIF
          ENDIF
 730    CONTINUE
!
 800  CONTINUE
!
      RETURN
!     end of subroutine SWOEXA
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOINA (XC, YC, AC2, ACLOC, KGRPNT, DEPXY, CROSS,EXCPT)  40.86 30.50
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. AUTHORS
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!     40.86: Nico Booij
!
!  1. Update
!
!     10.10, Aug. 94: separated from subr. SWOEXA
!     30.50,        : If depth on one of the corners is negative value 0
!                     is returned
!     30.72, Sept 97: Replaced DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     40.13, Aug. 01: provision for repeating grid (KREPTX>0)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.86, Feb. 08: modification to prevent interpolation over an obstacle
!
!  2. Purpose
!
!       interpolates local action density ACLOC from array AC2
!
!  4. Argument list
!
!       XC, YC  real   input    comp. grid coordinates
!       AC2     real a input    action densities
!       ACLOC   real a outp     local action density spectrum
!
!  5. SUBROUTINES CALLING
!
!       SWOEXA (SWAN/OUTP)
!
! 10. SOURCE TEXT
!
      LOGICAL :: EXCPT     ! if true value is undefined                   40.86
      LOGICAL :: CROSS(4)  ! true if obstacle is between output point     40.86
                           ! and computational grid point                 40.86
      REAL     XC, YC, AC2(MDC,MSC,MCGRD), ACLOC(MDC, MSC),               30.21
     &         DEPXY(MCGRD)
!
      INTEGER  KGRPNT(MXC,MYC)                                            30.21
!
      REAL :: WW(1:4)    ! Interpolation weights for the 4 corners        40.86
      REAL :: SUMWW      ! sum of the weights                             40.86
      INTEGER :: INDX(1:4)     ! grid counters for the 4 corners          40.86
      INTEGER :: JX(1:4), JY(1:4)  ! grid counters for the 4 corners      40.86
      INTEGER :: JC            ! corner counter                           40.86
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWOINA')
!
      EXCPT = .FALSE.                                                     40.86
      WW(1:4) = 0.                                                        40.86
      SUMWW = 0.                                                          40.86
      ACLOC = 0.                                                          40.86
      JX1 = INT(XC+3.) - 2
      JX2 = JX1+1
      SX2 = XC + 1. - FLOAT(JX1)
      SX1 = 1. - SX2
      IF (KREPTX .EQ. 0) THEN                                             40.13
        IF (SX1.LT.0.01 .OR. JX1.EQ.0)   THEN
          SX1 = 0.
          SX2 = 1.
          JX1 = MAX(1,JX1)
        ENDIF
        IF (SX2.LT.0.01 .OR. JX1.EQ.MXC) THEN
          SX2 = 0.
          SX1 = 1.
          JX2 = MIN(MXC,JX2)
        ENDIF
      ELSE                                                                40.13
!       repeating grid                                                    40.13
        JX1 = 1 + MODULO (JX1-1, MXC)                                     40.13
        JX2 = 1 + MODULO (JX2-1, MXC)                                     40.13
      ENDIF                                                               40.13
!
      IF (ONED) THEN                                                      40.86
        JY1 = 1                                                           40.86
        JY2 = 1                                                           40.86
        SY1 = 0.5                                                         40.86
        SY2 = 0.5                                                         40.86
      ELSE
        JY1 = INT(YC+3.) - 2
        JY2 = JY1+1
        SY2 = YC + 1. - FLOAT(JY1)
        SY1 = 1. - SY2
        IF (SY1.LT.0.01 .OR. JY1.EQ.0) THEN
          SY1 = 0.
          SY2 = 1.
          JY1 = MAX(1,JY1)
        ENDIF
        IF (SY2.LT.0.01 .OR. JY1.EQ.MYC) THEN
          SY2 = 0.
          SY1 = 1.
          JY2 = MIN(MYC,JY2)
        ENDIF
      ENDIF
!
!      *** Using indirect addressing for AC2   ***
!
      DO 91 ISIGM = 1, MSC                                                30.72
        DO 90 ID  = 1, MDC
          ACLOC(ID,ISIGM) = 0.
  90    CONTINUE                                                          30.72
  91  CONTINUE                                                            30.72
!
      IF (.NOT.EXCPT) THEN
         JX(1) = JX1                                                      40.86
         JY(1) = JY1                                                      40.86
         WW(1) = SX1*SY1                                                  40.86
         JX(2) = JX2                                                      40.86
         JY(2) = JY1                                                      40.86
         WW(2) = SX2*SY1                                                  40.86
         JX(3) = JX1                                                      40.86
         JY(3) = JY2                                                      40.86
         WW(3) = SX1*SY2                                                  40.86
         JX(4) = JX2                                                      40.86
         JY(4) = JY2                                                      40.86
         WW(4) = SX2*SY2                                                  40.86
         DO JC = 1, 4                                                     40.86
           INDX(JC) = KGRPNT(JX(JC),JY(JC))                               40.86 30.21
           IF (WW(JC).LT.0.01) THEN
             WW(JC) = 0.
           ELSE
             IF (INDX(JC).LE.1) THEN                                      40.86
               WW(JC) = 0.                                                40.86
             ELSE IF (DEPXY(INDX(JC)).LE.DEPMIN) THEN                     40.86
!              dry point
               EXCPT =  .TRUE.                                            40.94 40.86
             ELSE IF (CROSS(JC) .AND. WW(JC).LT.0.999) THEN               40.86
!              obstacle                                                   40.86
               WW(JC) = 0.                                                40.86
             ENDIF
           ENDIF
         ENDDO
         SUMWW = SUM(WW(1:4))                                             40.86
         IF (.NOT.EXCPT) THEN
           IF (SUMWW.GT.0.01) THEN                                        40.86
             DO JC = 1, 4
               IF (WW(JC).GT.1.E-6) THEN
                 DO ISIGM = 1, MSC
                   DO ID = 1, MDC
                     ACLOC(ID,ISIGM) = ACLOC(ID,ISIGM) +
     &                        WW(JC)*AC2(ID,ISIGM,INDX(JC))
                   ENDDO
                 ENDDO
               ENDIF
             ENDDO
             IF (SUMWW.LT.0.999999) THEN
               DO ISIGM = 1, MSC
                 DO ID = 1, MDC
                   ACLOC(ID,ISIGM) = ACLOC(ID,ISIGM) / SUMWW
                 ENDDO
               ENDDO
             ENDIF
           ELSE
             EXCPT =  .TRUE.
           ENDIF
         ENDIF
      ENDIF
      IF (ITEST.GE. 10) WRITE (PRTEST, 89)
     &   XC, YC, (JX(JC), JY(JC), WW(JC), INDX(JC), CROSS(JC), JC=1,4),
     &   SUMWW
  89  FORMAT (' SWOINA ', 2F9.3, 4(2X, 2I5, F6.3, 1X, I4, 1X, L1), 2X,
     &                    F6.3)
 900  RETURN
!     end of subroutine SWOINA
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOEXF (MIP      ,XC       ,YC       ,VOQR     ,
     &                   VOQ      ,AC2      ,DEP2     ,SPCSIG   ,         30.72
     &                   WK       ,CG       ,SPCDIR   ,NE       ,
     &                   NED      ,KGRPNT   ,XCGRID   ,YCGRID   ,         30.72
     &                   IONOD                                            40.31
     &                                                          )
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
!
      IMPLICIT NONE                                                       30.81
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.80, 40.13: Nico Booij
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.55, Mar. 97: Procedure updated for curvilinear coordinates basics is
!                     described in SWANDOC.WP5 comp. grid point coordinates are
!                     new arguments
!     30.72, Oct. 97: Logical function EQREAL introduced for floating point
!                     comparisons
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.80, Apr. 98: Provision for 1D computation
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Dec. 98: Argument list KSCIP1 adjusted
!     30.81, Dec. 98: Implicit none added, force for point surrounded by
!                     dry points set to 0.
!     40.13, Aug. 01: provision for repeating grid (KREPTX>0)
!                     spherical coordinates taken into account
!                     swcomm2.inc reactivated
!     40.31, Jan. 04: adapted for parallelisation with MPI
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculates wave-driven force (output quantity IVTYPE=20)
!
!  3. Method
!
!     Radiation stresses are defined as:
!                     /
!     Sxx = rho grav | ((N cos^2(theta) + N - 1/2) sig Ac) d sig d theta
!                   /
!                     /
!     Sxy = rho grav | (N sin(theta) cos(theta) sig Ac) d sig d theta
!                   /
!                     /
!     Syy = rho grav | ((N sin^2(theta) + N - 1/2) sig Ac) d sig d theta
!                   /
!
!     The force in x-direction and y-direction are:
!
!     Fx = - (@Sxx/@x + @Sxy/@y)
!     Fy = - (@Sxy/@x + @Syy/@y)
!
!     where @ denotes the partial derivative.
!
!     The value of N and its derivative w.rt. the depth are calculated
!     in the KSCIP1 subroutine.
!
!     First the gradients with respect to i and j (comp. grid counters)
!     are computed, then these are transformed into gradients in (x,y)
!
!  4. Argument variables
!
!     AC2     input  action density
!     CG      local  group velocity in output point
!     DEP2    input  depth at comp. grid points
!     KGRPNT  input  index for indirect adressing
!     IONOD   input  array indicating in which subdomain output           40.51
!                    points are located                                   40.31
!     MIP     input  number of output points
!     NE      local  ratio of group and phase velocity
!     NED     local  derivative of NE with respect to depth
!     SPCDIR  input  (*,1); spectral directions (radians)                 30.82
!                    (*,2); cosine of spectral directions                 30.82
!                    (*,3); sine of spectral directions                   30.82
!                    (*,4); cosine^2 of spectral directions               30.82
!                    (*,5); cosine*sine of spectral directions            30.82
!                    (*,6); sine^2 of spectral directions                 30.82
!     SPCSIG  input  relative frequencies in computational domain in
!                    sigma-space                                          30.72
!     XC, YC  input  comp. grid coordinates of output point
!     XCGRID  input  coordinates of computational grid in x-direction     30.72
!     YCGRID  input  coordinates of computational grid in y-direction     30.72
!     VOQR    input  location in VOQ of a certain outp quant.
!     VOQ     output values of output quantities
!     WK      local  wavenumber in output point
!
      INTEGER MIP, VOQR(*) ,KGRPNT(MXC,MYC)                               30.21
      INTEGER IONOD(*)                                                    40.31
      REAL    AC2(MDC,MSC,MCGRD), CG(*), DEP2(MCGRD), NE(*), NED(*)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XC(MIP), YC(MIP)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)                            30.72
      REAL    VOQ(MIP,*), WK(*)
      LOGICAL    EQREAL                                                   30.72
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     AC2LOC       Local action density
!     ACWAVE       Action density in output point
!     ACWI, ACWJ   dAC/dI and dAC/dJ
!     ACWX         X-gradient of local action density
!     ACWY         Y-gradient of local action density
!     DDET         determinant
!     DDI, DDJ     dDEP/dI and dDEP/dJ
!     DDX, DDY     spatial depth gradients
!     DIX, DIY     coefficients for transformation from I-gradient
!                  to (X,Y)-gradients
!     DJX, DJY     coefficients for transformation from J-gradient
!                  to (X,Y)-gradients
!     DS2
!     DXI, DXJ     dX/dI and dX/dJ
!     DYI, DYJ     dY/dI and dY/dJ
!     DEP          depth
!     DEPLOC       local depth
!     FX, FY       (preliminary) forces in X-, and Y-direction
!     FXADD, FYADD Cumulated forces/(RHO*GRAV) per frequency and directi-
!                  onal step in X-, and Y-direction
!     ID           counter for steps in direction
!     IENT         number of entries
!     IND1, IND2,
!     IND3, IND4,
!     IND5, IND6,
!     IND7, IND8,
!     IND9         indirect adresses
!     IP           counter
!     IS           counter for sigma
!     IVTYPE
!     JX           counter in X-direction
!     JXLO, JXUP   lower resp. upper gridpoint number of point under consideration
!                  in X-direction
!     JY           counter in Y-direction
!     JYLO, JYUP   lower resp. upper gridpoint number of point under consideration
!                  in Y-direction
!     NAX, NAY     derivative of N * Ac.dens. = N * E / Sigma, w.r.t. X
!                  or Y, respectively.
!     ONX, ONY     Indicates whether or not a output point lies on a
!                  computational point or not
!     RRDI,RRDJ    multiplication factor: 0.5 in case of two-sided or 1 in case
!                  of one-sided differential
!     SIG          dummy variable
!     SXLO, SXUP   weight coefficients for the lower and upper x-level of the
!                  point under consideration, respectively.
!     SYLO, SYUP   weight coefficients for the lower and upper y-level of the
!                  point under consideration, respectively.
!
      REAL        ACWAV, ACWI, ACWJ, ACWX, ACWY, DDET, DDI, DDJ, DDX,
     &            DDY, DIX, DIY,DJX, DJY, DS2, DXI, DXJ, DYI, DYJ, DEP,
     &            DEPLOC, FX,FY, FXADD, FYADD, NAX, NAY, RRDI, RRDJ,
     &            SIG, SXLO, SXUP, SYLO, SYUP, CSLAT
      INTEGER     ID, IENT, IND1, IND2, IND3, IND4, IND5, IND6, IND7,
     &            IND8, IND9, IP, IS, IVTYPE, JX, JXLO, JXUP, JY,
     &            JYLO, JYUP
      LOGICAL     ONX, ONY
      REAL        AC2LOC(MCGRD)                                           40.31
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     KSCIP1           calculates WK, CG, N and ND
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWEXCHG          exchanges AC2 at subdomain boundaries
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     -In determining derivatives one-sided differences are used at
!      border meshes; for output points inside a mesh derivative
!      over one step is taken; for output points on a computational
!      grid point a central derivative is taken
!     -A marigin of 0.01 m is taken outside the computational grid.
!     -The range of the counter runs from 1 to MXC; the range of XC(IP) runs
!      from 0 to MXC-1!
!
!  Counter: 1          2          JX        JX+1      JX+2       MXC-1       MXC
!
!         |=|----------|-- -- -- -|--------|=|=|--------|-- -- -- -|----------|=|
!
!  XC:      0          1                     JX                  MXC-2      MXC-1
!
!
!     -Order in which they are treated:
!
!         |=|----------|          |--------|=|=|--------|          |----------|=|
! Order:         A                     B     C      D                    E
!
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     For all output points do
!         Initialize both force components as 0
!         Determine neighbouring points to be used for gradients in
!         X and Y
!         Call KSCIP1 (determine derivative of N with respect to depth)
!         For all spectral components do
!             determine derivative of nE with respect to X
!             determine derivative of nE with respect to Y
!             calculate contribution to force components
!     ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'SWOEXF')
!
      IVTYPE = 20
!
!     --- exchange action densities at subdomain interfaces               40.31
!         within distributed-memory environment                           40.31
!
      IF (PARLL) THEN                                                     40.41
!TIMG         CALL SWTSTA(213)                                                 40.31
         DO ID = 1, MDC                                                   40.31
            DO IS = 1, MSC                                                40.31
               AC2LOC(:) = AC2(ID,IS,:)                                   40.31
               CALL SWEXCHG( AC2LOC, KGRPNT )                             40.31
               AC2(ID,IS,:) = AC2LOC(:)                                   40.31
            END DO                                                        40.31
         END DO                                                           40.31
!TIMG         CALL SWTSTO(213)                                                 40.31
         IF (STPNOW()) RETURN                                             40.31
      END IF
!
!     loop over all output points
!
      DO 800 IP=1,MIP
        DEP = VOQ(IP,VOQR(4))
        IF (DEP.LE.0.)                  GOTO 700
        IF (EQREAL(DEP,OVEXCV(4)))      GOTO 700                          30.72
        IF (KREPTX .EQ. 0) THEN                                           40.13
          IF (XC(IP).LT.-0.01)            GOTO 700
          IF (XC(IP).GT.REAL(MXC-1)+0.01) GOTO 700
        ENDIF                                                             40.13
        IF (YC(IP).LT.-0.01)            GOTO 700
        IF (YC(IP).GT.REAL(MYC-1)+0.01) GOTO 700
        IF (PARLL .AND. IONOD(IP).NE.INODE) GOTO 700                      40.31
!
!       first the action density spectrum is interpolated
!
        FX  = 0.
        FY  = 0.
        JX  = NINT(XC(IP))
        RRDI = 1.
        ONX = .FALSE.
        IF (KREPTX.EQ.0 .AND. JX.EQ.0) THEN                               40.13
          JXLO = 1
          JXUP = 2
          JX   = 1
          SXUP = XC(IP)
          SXLO = 1.-SXUP
        ELSE IF (KREPTX.EQ.0 .AND. JX.EQ.MXC-1) THEN                      40.13
          JXLO = MXC-1
          JXUP = MXC
          SXLO = REAL(MXC-1)-XC(IP)
          SXUP = 1.-SXLO
          JX   = JX+1                                                     30.81
        ELSE IF (XC(IP).LT.REAL(JX)-0.01) THEN
          JXLO = JX
          JXUP = JX+1
          SXLO = REAL(JX)-XC(IP)
          SXUP = 1.-SXLO
          JX   = JX+1                                                     30.81
        ELSE IF (XC(IP).GT.REAL(JX)+0.01) THEN
          JXLO = JX+1
          JXUP = JX+2
          SXUP = XC(IP)-REAL(JX)
          SXLO = 1.-SXUP
          JX   = JX+1                                                     30.81
        ELSE
          JXLO = JX
          JXUP = JX+2
          RRDI = 0.5
          JX   = JX+1
          ONX  = .TRUE.
        ENDIF
        IF (KREPTX .GT. 0) THEN                                           40.13
          JX   = 1 + MODULO (JX-1, MXC)                                   40.13
          JXLO = 1 + MODULO (JXLO-1, MXC)                                 40.13
          JXUP = 1 + MODULO (JXUP-1, MXC)                                 40.13
        ENDIF                                                             40.13
        IF (ONED) THEN                                                    30.80
          JYLO = 1                                                        30.80
          JYUP = 1                                                        30.80
          JY   = 1                                                        30.80
          RRDJ = 0.                                                       30.80
          ONY  = .TRUE.                                                   30.80
        ELSE                                                              30.80
          JY   = NINT(YC(IP))
          RRDJ = 1.
          ONY  = .FALSE.
          IF (JY.EQ.0) THEN
            JYLO = 1
            JYUP = 2
            JY   = 1
            SYUP = YC(IP)
            SYLO = 1.-SYUP
          ELSE IF (JY.EQ.MYC-1) THEN
            JYLO = MYC-1
            JYUP = MYC
            SYLO = REAL(MYC-1)-YC(IP)
            SYUP = 1.-SYLO
            JY   = JY+1                                                   30.81
          ELSE IF (YC(IP).LT.REAL(JY)-0.01) THEN
            JYLO = JY
            JYUP = JY+1
            SYLO = REAL(JY)-YC(IP)
            SYUP = 1.-SYLO
            JY   = JY+1                                                   30.81
          ELSE IF (YC(IP).GT.REAL(JY)+0.01) THEN
            JYLO = JY+1
            JYUP = JY+2
            SYUP = YC(IP)-REAL(JY)
            SYLO = 1.-SYUP
            JY   = JY+1                                                   30.81
          ELSE
            JYLO = JY
            JYUP = JY+2
            RRDJ = 0.5
            JY  = JY+1
            ONY = .TRUE.
          ENDIF
        ENDIF                                                             30.80
!
!       *** Using indirect addressing for arrays AC2 and DEP2 ***
        IND1 = KGRPNT(JXLO,JYLO)                                          30.21
        IND2 = KGRPNT(JXUP,JYLO)                                          30.21
        IND3 = KGRPNT(JXUP,JYUP)                                          30.21
        IND4 = KGRPNT(JXLO,JYUP)                                          30.21
        IND5 = KGRPNT(JXLO,JY  )                                          30.21
        IND6 = KGRPNT(JXUP,JY  )                                          30.21
        IND7 = KGRPNT(JX  ,JYLO)                                          30.21
        IND8 = KGRPNT(JX  ,JYUP)                                          30.21
        IND9 = KGRPNT(JX  ,JY  )                                          30.21
        IF (ONY) THEN                                                     40.00
          IF (DEP2(IND5).LE.DEPMIN) GOTO 700
          IF (DEP2(IND6).LE.DEPMIN) GOTO 700
        ELSE
          IF (DEP2(IND1).LE.DEPMIN) GOTO 700
          IF (DEP2(IND2).LE.DEPMIN) GOTO 700
          IF (DEP2(IND3).LE.DEPMIN) GOTO 700
          IF (DEP2(IND4).LE.DEPMIN) GOTO 700
        ENDIF
        IF (ONX) THEN                                                     40.00
          IF (DEP2(IND7).LE.DEPMIN) GOTO 700
          IF (DEP2(IND8).LE.DEPMIN) GOTO 700
        ELSE
          IF (DEP2(IND1).LE.DEPMIN) GOTO 700
          IF (DEP2(IND2).LE.DEPMIN) GOTO 700
          IF (DEP2(IND3).LE.DEPMIN) GOTO 700
          IF (DEP2(IND4).LE.DEPMIN) GOTO 700
        ENDIF
!
!       determine depth and (x,y) derivatives w.r.t. i and j
!
        IF (ONY) THEN
          DDI = RRDI * (DEP2(IND6)-DEP2(IND5))
          DXI = RRDI * (XCGRID(JXUP,JY)-XCGRID(JXLO,JY))                  30.72
          DYI = RRDI * (YCGRID(JXUP,JY)-YCGRID(JXLO,JY))                  30.72
        ELSE
          DDI = RRDI * (SYUP*(DEP2(IND3)-DEP2(IND4)) +
     &                  SYLO*(DEP2(IND2)-DEP2(IND1)))
          DXI = RRDI * (SYUP*(XCGRID(JXUP,JYUP)-XCGRID(JXLO,JYUP)) +      30.72
     &                  SYLO*(XCGRID(JXUP,JYLO)-XCGRID(JXLO,JYLO)))       30.72
          DYI = RRDI * (SYUP*(YCGRID(JXUP,JYUP)-YCGRID(JXLO,JYUP)) +      30.72
     &                  SYLO*(YCGRID(JXUP,JYLO)-YCGRID(JXLO,JYLO)))       30.72
        ENDIF
        IF (ONX) THEN
          DDJ = RRDJ * (DEP2(IND8)-DEP2(IND7))
          DXJ = RRDJ * (XCGRID(JX,JYUP)-XCGRID(JX,JYLO))                  30.72
          DYJ = RRDJ * (YCGRID(JX,JYUP)-YCGRID(JX,JYLO))                  30.72
        ELSE
          DDJ = RRDJ * (SXUP*(DEP2(IND3)-DEP2(IND2)) +
     &                  SXLO*(DEP2(IND4)-DEP2(IND1)))
          DXJ = RRDJ * (SXUP*(XCGRID(JXUP,JYUP)-XCGRID(JXUP,JYLO)) +      30.72
     &                  SXLO*(XCGRID(JXLO,JYUP)-XCGRID(JXLO,JYLO)))       30.72
          DYJ = RRDJ * (SXUP*(YCGRID(JXUP,JYUP)-YCGRID(JXUP,JYLO)) +      30.72
     &                  SXLO*(YCGRID(JXLO,JYUP)-YCGRID(JXLO,JYLO)))       30.72
        ENDIF
        IF (KSPHER.GT.0) THEN                                             40.13
!         spherical coordinates are used; first compute cos(latitude)     40.13
          CSLAT = COS(DEGRAD*(YOFFS+YCGRID(JX,JY)))                       40.61 40.13
!         LENDEG is the length of one degree of the sphere                40.13
          DXI = DXI * LENDEG * CSLAT                                      40.13
          DYI = DYI * LENDEG                                              40.13
          DXJ = DXJ * LENDEG * CSLAT                                      40.13
          DYJ = DYJ * LENDEG                                              40.13
        ENDIF                                                             40.13
!
!       coefficients from transformation from (i,j)-gradients to (x,y)-gradients
!
        IF (JXUP.EQ.JXLO .AND. JYUP.EQ.JYLO) THEN                         30.81
!         point surrounded by dry points                                  30.81
          DIX  = 0.                                                       30.81
          DIY  = 0.                                                       30.81
          DJX  = 0.                                                       30.81
          DJY  = 0.                                                       30.81
        ELSE IF (JXUP.EQ.JXLO) THEN                                       30.80
!         no forces in i-direction                                        30.81
          DS2  = DXJ**2 + DYJ**2                                          30.80
          DIX  = 0.                                                       30.80
          DIY  = 0.                                                       30.80
          DJX  = DXJ/DS2                                                  30.80
          DJY  = DYJ/DS2                                                  30.80
        ELSE IF (JYUP.EQ.JYLO) THEN                                       30.80
!         no forces in j-direction                                        30.81
          DS2  = DXI**2 + DYI**2                                          30.80
          DIX  = DXI/DS2                                                  30.80
          DIY  = DYI/DS2                                                  30.80
          DJX  = 0.                                                       30.80
          DJY  = 0.                                                       30.80
        ELSE                                                              30.80
!         coefficients for transformation from                            30.81
!         (i,j)-gradients to (x,y)-gradients                              30.81
          DDET = DXI*DYJ - DXJ*DYI
          DIX  =  DYJ / DDET
          DIY  = -DXJ / DDET
          DJX  = -DYI / DDET
          DJY  =  DXI / DDET
        ENDIF                                                             30.80
!       spatial depth gradients:
        DDX  = DDI*DIX + DDJ*DJX
        DDY  = DDI*DIY + DDJ*DJY
!
        IF (ITEST.GE.80 .OR. IOUTES .GE. 20) WRITE (PRTEST, 88) IP,
     &  JXLO, JXUP, JYLO, JYUP ,SXLO, SXUP, SYLO, SYUP,
     &  DIX, DIY, DJX, DJY                                                30.80
  88    FORMAT (' SWOEXF ', 5I6, 2X, 4F7.4, 2X, 4E12.4)                   30.80
!
!       compute NE and NED
!
        DEPLOC = VOQ(IP,VOQR(4))
        CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)                30.81 30.72
        IF (ITEST.GE.100 .OR. IOUTES .GE. 20) THEN
          WRITE (PRTEST, 98)  DEPLOC, DDX, DDY
  98      FORMAT (' depth & gradient ', 4(1X,F9.4))
          DO 100 IS = 1, MIN(MSC,20)
             WRITE (PRTEST, 99) IS, SPCSIG(IS), NE(IS),                   30.72
     &                          NED(IS)
  99         FORMAT (' i, SPCSIG, N, Nd ', I2, 3(1X, E12.4))              30.72
 100      CONTINUE
        ENDIF
!
        DO 300 ID  = 1, MDC
          DO 290 IS = 1, MSC                                              30.81 18/MAR
            SIG = SPCSIG(IS)                                              30.72
!
!           ACWAV is local action density
!
            IF (ONX.AND.ONY) THEN
               ACWAV = AC2(ID,IS,IND9)
            ELSE IF (ONX) THEN
               ACWAV = SYLO * AC2(ID,IS,IND7) +
     &                 SYUP * AC2(ID,IS,IND8)
            ELSE IF (ONY) THEN
               ACWAV = SXLO * AC2(ID,IS,IND5) +
     &                 SXUP * AC2(ID,IS,IND6)
            ELSE
               ACWAV = SXLO * (SYLO * AC2(ID,IS,IND1) +
     &                         SYUP * AC2(ID,IS,IND4)) +
     &                 SXUP * (SYLO * AC2(ID,IS,IND2) +
     &                         SYUP * AC2(ID,IS,IND3))
            ENDIF
!
!           ACWX is X-gradient of local action density, ACWY is Y-gradient
!
            IF (ONY) THEN
               ACWI = RRDI * (AC2(ID,IS,IND6) -
     &                        AC2(ID,IS,IND5))
            ELSE
               ACWI = RRDI * (SYLO * (AC2(ID,IS,IND2) -
     &                                AC2(ID,IS,IND1)) +
     &                        SYUP * (AC2(ID,IS,IND3) -
     &                                AC2(ID,IS,IND4)))
            ENDIF
            IF (ONX) THEN
               ACWJ = RRDJ * (AC2(ID,IS,IND8) -
     &                        AC2(ID,IS,IND7))
            ELSE
               ACWJ = RRDJ * (SXLO * (AC2(ID,IS,IND4) -
     &                                AC2(ID,IS,IND1)) +
     &                        SXUP * (AC2(ID,IS,IND3) -
     &                                AC2(ID,IS,IND2)))
            ENDIF
!
!           spatial action density gradients:                             30.55
            ACWX = ACWI*DIX + ACWJ*DJX
            ACWY = ACWI*DIY + ACWJ*DJY
!
!           NAX is the derivative of N * Ac.dens.  w.r.t. X
!           So NAX = @(N*Ac)/@X =Ac*@N/@X +N*@Ac/@X
!
!           where @ denotes the partial derivative.
!
!           Further note that that @N/@X = @N/@Depth * @Depth/@X
!                                        = NED * DDX
!
!           Anologously for NAY.
!
            NAX = NE(IS) * ACWX + NED(IS) * DDX * ACWAV
            NAY = NE(IS) * ACWY + NED(IS) * DDY * ACWAV
            FXADD = - ( (SPCDIR(ID,4) + 1.) * NAX - 0.5 * ACWX +          20.44
     &                   SPCDIR(ID,5) * NAY ) * SIG
            FYADD = - ( (SPCDIR(ID,6) + 1.) * NAY - 0.5 * ACWY +          20.44
     &                   SPCDIR(ID,5) * NAX ) * SIG
!
!           integration
!
            FX = FX + SIG * FXADD                                         20.35
            FY = FY + SIG * FYADD                                         20.35
 290      CONTINUE
 300    CONTINUE
!
        FX = RHO * GRAV * FX * DDIR * FRINTF                              20.77
        FY = RHO * GRAV * FY * DDIR * FRINTF                              20.77
        VOQ(IP,VOQR(IVTYPE))   = (COSCQ*FX - SINCQ*FY)
        VOQ(IP,VOQR(IVTYPE)+1) = (SINCQ*FX + COSCQ*FY)
        GOTO 800
!
!       points on land: assign exception value
!
 700    VOQ(IP,VOQR(IVTYPE))   = OVEXCV(IVTYPE)
        VOQ(IP,VOQR(IVTYPE)+1) = OVEXCV(IVTYPE)
!
 800  CONTINUE
!
      RETURN
!     end of subroutine SWOEXF
      END
