!
!     SWAN/READ file 1 of 2
!
! Contents of this file:
!
!     SWREAD:  Reading and processing of the user commands describing the model
!     SINPGR:  Read parameters of an input grid
!     SREDEP
!     SSFILL
!     CGINIT
!     SWDIM
!     CGBOUN   determines boundary of true computational region
!     INITVA:  Processing command INIT and compute initial state of       30.70
!              the wave field                                             30.70
!     BACKUP                                                              40.00
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWREAD (COMPUT)                                          40.31 30.90
!                                                                      *
!***********************************************************************

!     Modules

      USE TIMECOMM                                                        40.41
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.13
      USE M_SNL4                                                          40.17
      USE M_GENARR                                                        40.31
      USE M_OBSTA                                                         40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
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
!     30.60: Nico Booij
!     30.61: Roberto Padilla
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.75: Nico Booij
!     30.80: Nico Booij
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma
!     32.01: Roeland Ris & Cor van der Schelde
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     32.06: Roeland Ris
!     33.08: W. Erick Rogers
!     33.09: Nico Booij
!     33.10: W. Erick Rogers and Nico Booij
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.09: Annette Kieftenburg
!     40.10: IJsbrand Haagsma
!     40.12: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.14: Annette Kieftenburg
!     40.16: IJsbrand Haagsma
!     40.17: IJsbrand Haagsma
!     40.18: Annette Kieftenburg
!     40.21: Agnieszka Herman
!     40.23: Marcel Zijlema
!     40.28: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.38: Annette Kieftenburg
!     40.08: W. Erick Rogers
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     40.55: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.87: Marcel Zijlema
!
!  1. Updates
!
!     30.60, July 97: command CGRID, exception values added
!     30.60, Aug. 97: JCOOX and JCOOY used when addingco ordinate arrays change in
!                     command BOUND option SEGM
!     30.60, Aug. 97: PNUMS(20) is set to 0.1 in command WIND if third generation
!                     is used
!     30.60, Aug. 97: argument XYTST added in call of SWDIM
!     30.60, Aug. 97: activate initial condition in commands WIND and GEN3
!     30.60, Aug. 97: command CGRID, default keyword is REG keyword EXC is not
!                     required
!     30.60, Aug. 97: command OBST, names changed into ALPHA and BETA; control
!                     strings changed from UNC into STA
!     30.60, Aug. 97: uncommented statement CALL WRNKEY
!     30.70, Sep. 97: value of PWTAIL(1) is set to 5 in command GEN3 JANS,
!                     GROWTH JANS and WCAP JANS
!     30.72, Oct. 97: logical function EQREAL introduced for floating point
!                     comparisons
!     30.72, Nov. 97: Added the command syntax for all the commands as comments
!                     from the user manual
!     30.72, Nov. 97: Header renewed, updated method and argument variable
!                     description
!     30.72, Nov. 97: Did set the correct pointers for command GEN3 JANS, as
!                     was already done correctly for command WCAP JANS
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Jan. 98: Removed reference to quadruplets in WIND command
!     32.01, Jan. 98: Introduced SET NAUT command (project h3268)
!     32.01, Jan. 98: Introduced keyword CON/VAR in command
!                     BOU STAT ... SPEC1D/SPEC2D
!     30.72, Jan. 98: Moved BOU NONSTAT non operational warning to end of
!                     command BOU
!     30.72, Jan. 98: Changed size of JAUX(7) array to MCGRD
!     32.01, Jan. 98: Modifications for nautical convention, interpolation of
!                     spectra at the boundary and warning (project h3268)
!     32.02, Feb. 98: Introduced 1D-version
!     30.70, Feb. 98: MODE DYN modified into MODE NONSTationary
!                     option WINDGrowth added in command OFF
!                     Nautical convention introduced into command CGRID (sector)
!                     in command SET name 'negmes' changed into 'maxmes'
!     30.72, Mar. 98: Leave limiter on when ITRIAD > 0 in command OFF QUAD
!     30.70, Mar. 98: option CON/VAR after options SPEC1D/SPEC2D
!                     keyword STAT made optional
!                     name GRWMX changed into LIMITER
!                     assignment of limiter in command OFF QUAD corrected
!     30.75, Mar. 98: set ICOND=1 (default init.cond.) for nonstationary mode
!     30.82, Apr. 98: removed reference to commons KAART and KAR
!     40.00, Sep. 98: in command OBST square of TRCOEF is stored (this is used
!                     as transmission coefficient for action density)
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.81, Nov. 98: Adjustment for 1-D case of new boundary conditions
!     30.80, Nov. 98: Provision for limitation on Ctheta (refraction)
!     40.00, Jan. 99: new command OUTPUT QUANTITY: allows user to change
!                     properties of output quantities.
!     34.01, Feb. 99: Introducing STPNOW
!     30.82, Mar. 99: Deactivate limiter for GEN1 and GEN2
!     40.00, Apr. 99: restructure command MODE: Nonstat Oned is now possible
!     33.08, July 98: input for model with higher order "S&L" scheme.
!     33.09, Aug. 98: input for spherical coordinates
!     32.06, June 99: Set correct values for IGEN
!     30.82, Aug. 99: Modified default values for PTRIAD, after deactivating
!                     limiter for only the triads.
!     30.82, Aug. 99: Modified command NUM to include settings for the SETUP
!                     and the global stop criterion
!     40.01, Sep. 99: XASM and YASM replace fixed numbers
!     40.03, Dec. 99: command QUANTITY corrected
!     40.10, Mar. 00: prepared for exact quadruplets
!     33.10, Jan. 00: input for model with higher order "SORDUP" scheme
!     40.03, May  00: command INCLude added, array INCNUM added
!     40.09, May  00: TRCOEF**2 replaced by TRCOEF (to make command options
!                     TRANSM and DAM consistent with one another in subroutine
!                     SWTRCOEF in swanser)
!     40.09, May. 00: Reflection command option added
!     40.03, Aug. 00: error message if start time is before current time (command COMP)
!            Sep. 00: inconsistency with manual corrected
!     40.02, Sep. 00: Changed SCHEME command to PROP and handling [cdlim] modified
!     40.02, Oct. 00: In case of Nautical directions then output in not w.r.t. frame
!     40.02, Oct. 00: Recalculate whitecapping coefficients for new SWCAP routine
!     40.02, Oct. 00: Avoided real/int conflict by introducing RPOOL in SREDEP
!     40.02, Oct. 00: Avoided real/int conflict by introducing replacing
!                     RPOOL for POOL in various calls
!     40.02, Oct. 00: Initialisation of IERR
!     40.12, Feb. 01: Avoided type conflict for OUTPS
!     40.18, Apr. 01: Reflection option extended
!     40.13, July 01: reading of PTRIAD(4) added in command TRIad
!                     command OUTPut OPTions added; module OUTP_DATA added
!                     command OUTPUT QUANTITY is made obsolete
!     40.13, Aug. 01: [xpc] and [ypc] are required in case of spherical coordinates
!                     [ylenc] and [myc] not required in 1-D mode
!     40.13, Oct. 01: USE OUTP_DATA added in view of longer filenames
!     40.13, Oct. 01: value of SPDIR1 in full circle case changed
!     40.13, Nov. 01: size of array BSPAUX increased
!     40.13, Nov. 01: command OUTPUT OPTIONS added
!     40.18, Apr. 01: Reflection option extended, scatter
!     40.28, Dec. 01: Reflection option extended, freq. dep.
!     40.38, Feb. 02: Reflection option extended, diffuse
!     40.14, Dec. 01: Extra check on user defined coefficients added
!     40.16, Dec. 01: Implemented limiter switches
!     40.17, Dec. 01: Implemented Multiple DIA
!     40.21, Aug. 01: Diffraction option added
!     40.23, Aug. 02: under-relaxation factor added
!     40.23, Sep. 02: coefficient urslim must be stored in PTRIAD(5)
!     40.23, Nov. 02: keyword FLUXLIM added
!     40.30, Feb. 03: introduction distributed-memory approach using MPI
!     40.08, Mar. 03: warning message added for use of curvilinear
!                     coordinate system
!     40.31, Dec. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Jun. 06: correct input of MDIA
!     40.55, Dec. 05: introducing vegetation model
!     40.80, Jun. 07: extension to unstructured grids
!     40.87, Apr. 08: in keyword QUANTITY more than 1 output parameter at once
!                     is possible and addition of [fmin] and [fmax]
!
!  2. Purpose
!
!     Reading and processing of the user commands describing the model
!
!  3. Method (updated 30.72)
!
!     A new line is read, in which the first keyword determines what the
!     command is. The command is read and processed. Common variable are
!     given proper values. After processing
!     the command the program returns to label 100, to process a new command.
!     This is repeated until the command STOP is found or the end of file is
!     reached.
!
!     Depending on the commands, the argument variable COMPUT is given a value,
!     depending on which the program will make a computation is some form or not.
!
!  4. Argument variables (updated 30.72)
!
!     COMPUT   Output variable that determine the sort of computation to be
!              performed by SWAN
!              ='COMP'; computation requested
!              ='NOCO'; no computation but output requested
!              ='RETR'; retrieve data previous computation
!              ='STOP'; make computation, output and stop
!
!  5. Parameter variables
!
      INTEGER, PARAMETER :: MXINCL = 10                                   40.03
      INTEGER, PARAMETER :: NVOTP  = 15                                   40.87
!
!  6. Local variables
!
!     DUM       Dummy variable                                            40.18
!     FD1       Coeff. for freq. dep. reflection: vertical displacement   40.28
!     FD2       Coeff. for freq. dep. reflection: shape parameter         40.28
!     FD3       Coeff. for freq. dep. reflection: directional coefficient 40.28
!     FD4       Coeff. for freq. dep. reflection: bending point of freq.  40.28
!     ICNL4     Counter for reading CNL4_? values                         40.17
!     ILAMBDA   Counter for reading LAMBDA values                         40.17
!     INCLEV    Include level, increases at INCL command,                 40.03
!               decreases at end-of-file                                  40.03
!     INCNUM    Unit reference numbers of included files                  40.03
!     ITMP1     auxiliary integer                                         40.31
!     ITMP2     auxiliary integer                                         40.31
!     ITMP3     auxiliary integer                                         40.31
!     ITMP4     auxiliary integer                                         40.31
!     LREF      Indicates whether reflection is active (#0.) or not (=0.) 40.09
!     LREFDIFF  Indicates whether scattered reflection is active (#0.)    40.18
!               or not (=0.)                                              40.18
!     LRFRD     Indicates whether frequency dependent reflection is       40.28
!               active (#0.) or not (=0.)                                 40.28
!     MORE      Indicates whether more entries of a loop need to be       40.17
!               performed                                                 40.17
!     POWN      user defined power of redistribution function             40.18
!     RLAMBDA   Dummy array containing lambda values for quadruplets      40.17
!
      INTEGER, SAVE     :: INCNUM(1:MXINCL) = 0                           40.03
      INTEGER, SAVE     :: INCLEV = 1                                     40.03
!
      INTEGER           :: IOSTAT = 0                                     40.02
      INTEGER           :: ICNL4, ILAMBDA                                 40.17
      INTEGER           :: ITMP1, ITMP2, ITMP3, ITMP4                     40.31
!
      LOGICAL           :: MORE                                           40.17
      LOGICAL, SAVE     :: LOBST = .FALSE.                                40.31
!
      INTEGER        :: LREF, LREFDIFF, LRFRD                             40.31
      REAL           :: POWN, DUM                                         40.31 40.18
      REAL           :: FD1, FD2, FD3, FD4                                40.31 40.28
      REAL, ALLOCATABLE :: RLAMBDA(:)                                     40.17

      CHARACTER*6  QOVSNM                                                 40.87
      CHARACTER*40 QOVLNM                                                 40.87
      REAL         QR(9)                                                  40.87

      TYPE(OBSTDAT), POINTER :: OBSTMP                                    40.31
      TYPE(OBSTDAT), SAVE, POINTER :: COBST                               40.31

      TYPE(OPSDAT), POINTER :: OPSTMP                                     40.31

      TYPE XYPT                                                           40.31
        REAL                :: X, Y
        TYPE(XYPT), POINTER :: NEXTXY
      END TYPE XYPT

      TYPE(XYPT), TARGET  :: FRST                                         40.31
      TYPE(XYPT), POINTER :: CURR, TMP                                    40.31

      TYPE AUXT                                                           40.87
        INTEGER             :: I                                          40.87
        TYPE(AUXT), POINTER :: NEXTI                                      40.87
      END TYPE AUXT                                                       40.87
      TYPE(AUXT), TARGET  :: FRSTQ                                        40.87
      TYPE(AUXT), POINTER :: CURRQ, TMPQ                                  40.87

      TYPE VEGPT                                                          40.55
        INTEGER              :: N                                         40.55
        REAL                 :: H, D, C                                   40.55
        TYPE(VEGPT), POINTER :: NEXTV                                     40.55
      END TYPE VEGPT                                                      40.55

      TYPE(VEGPT), TARGET  :: FRSTV                                       40.55
      TYPE(VEGPT), POINTER :: CURRV, TMPV                                 40.55
!
!  8. Subroutines used
!
!     DEGCNV: Transforms dir. from nautical to cartesian or vice versa    32.01
!     INITVA: Processing comm. INIT and comp. initial state of wave field 30.70
!     SINPGR: Read parameters of an input grid
!     SWINIT
!     REINC
!     SWRBC
!     SREDEP
!     SPRCON
!     SPROUT
!     RETSTP  read test points                                            40.00
!     SWNDPR (SWAN/SWREAD)
!     OCPINI
!     NWLINE
!     INKEYW
!     INCSTR
!     ININTG
!     INREAL
!     FOR
!     KEYWIS
      LOGICAL :: KEYWIS                                                   40.03
!     COPYCH
      LOGICAL :: EQREAL                                                   40.00
!     (all Ocean Pack)
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     The description of the structure of this subroutine is very
!     short as most of the source code can easily be understood with
!     the aid of the command descriptions in the user manual and the
!     purpose of the subroutines from the system documentation.
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     Call NWLINE for reading new line of user input
!     Call INKEYW to read a new command from user input
!     If the command is equal to one of the SWAN commands, then
!         Read and process the rest of the command
!     ----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL :: FOUND                                                    40.03
      LOGICAL, SAVE :: RUNMADE = .FALSE.                                  40.03
!
!     *** The logical variable LOGCOM has a record about which    ***
!     *** commands have been given to know if all the information ***
!     *** for certain command is available:                       ***     40.31
!     *** LOGCOM(2): command CGRID has been carried out           ***     40.31
!     *** LOGCOM(3): command READINP BOTTOM has been carried out  ***     40.31
!     *** LOGCOM(4): command READ COOR has been carried out       ***     40.31
!     *** LOGCOM(5): command READ UNSTRUC has been carried out    ***     40.80
!     *** LOGCOM(6): array AC2 has been allocated                 ***     40.31
!     *** In the current version, LOGCOM(1) has no meanings       ***     40.80 40.31
!
      LOGICAL, SAVE :: LOGCOM(1:6) = .FALSE.                              40.03
      INTEGER   TIMARR(6)                                                 40.31 30.00
      CHARACTER PSNAME *8, PNAME *8, COMPUT *(*), PTYPE *1, DTTIWR *18    40.00
      INTEGER, SAVE :: IENT = 0       ! number of entries to this subr
      INTEGER, SAVE :: LWINDR = 0     ! if non-zero, there is wind
      INTEGER, SAVE :: LWINDM = 3     ! type of wind growth formulation
      INTEGER, ALLOCATABLE :: IARR(:)                                     40.31
      INTEGER IVOTP(NVOTP), INDX(1)                                       40.87
      DATA IVOTP /10, 11, 13, 15, 16, 17, 18, 19,                         40.87
     &            28, 32, 33, 42, 43, 47, 48    /                         40.87

      CALL STRACE (IENT, 'SWREAD')
!
!     ***** read command *****
!
 100  CALL NWLINE
      IF (ELTYPE.EQ.'EOF') THEN
!       end-of-file encountered in (included) input file                  40.03
!       return to previous input file
        CLOSE (INCNUM(INCLEV))                                            40.03
        INCLEV = INCLEV - 1
        IF (INCLEV.EQ.0) THEN
          CALL MSGERR (4, ' unexpected end of command input')             40.03
          RETURN
        ENDIF
        ELTYPE = 'USED'
        INPUTF = INCNUM(INCLEV)                                           40.03
      ENDIF
!
      IF ( ITEST .GE. 200) THEN                                           30.70
        WRITE (PRTEST,*) ' BNAUT NA LABEL 100 IN SWREAD =', BNAUT         32.01
      ENDIF                                                               32.01
!
      CALL INKEYW ('REQ',' ')
!
!     ------------------------------------------------------------------
!                 PROCESSING OF COMMANDS
!     ------------------------------------------------------------------
!
!     STOP                                                                30.21
!
! ============================================================
!
! STOP
!
! ============================================================
!
      IF (KEYWIS ('STOP')) THEN
        IF (RUNMADE) THEN
          COMPUT = 'STOP'
        ELSE
          WRITE (PRINTF,*)' ** No computation requested **'
          COMPUT = 'NOCO'
        ENDIF
        RETURN
      ENDIF
!
!     ------------------------------------------------------------------
!
!     PROJECT      reading of project title and description
!
! ===============================================================
!
! PROJect  'NAME'  'NR'
!
!          'title1'
!
!          'title2'
!
!          'title3'
!
! ===============================================================
!
      IF (KEYWIS ('PROJ')) THEN
        CALL INCSTR ('NAME', PROJID, 'UNC', BLANK)
        CALL INCSTR ('NR', PROJNR, 'REQ', BLANK)
        CALL NWLINE
        IF (STPNOW()) RETURN                                              34.01
        CALL INCSTR ('TITLE1',PROJT1,'UNC',' ')
        CALL NWLINE
        IF (STPNOW()) RETURN                                              34.01
        CALL INCSTR ('TITLE2',PROJT2,'UNC',' ')
        CALL NWLINE
        IF (STPNOW()) RETURN                                              34.01
        CALL INCSTR ('TITLE3',PROJT3,'UNC',' ')
        GOTO 100
      ENDIF
!
!     INCLUDE     include another file in command file
!     ------------------------------------------------------------------
!     INCLude  'FILE'
!     ------------------------------------------------------------------
!
      IF (KEYWIS ('INCL')) THEN
        IF (INCLEV.EQ.1) INCNUM(INCLEV) = INPUTF                          40.03
        CALL INCSTR ('FILE' , FILENM , 'REQ', ' ')
        INCLEV = INCLEV + 1
        IF (INCLEV.GT.MXINCL) THEN
          CALL MSGERR (4, 'too many INCLUDE levels')                      40.03
          RETURN
        ENDIF
        IOSTAT = 0
        CALL FOR (INCNUM(INCLEV), FILENM, 'OF', IOSTAT)                   40.03
        INPUTF = INCNUM(INCLEV)
        GOTO 100
      ENDIF
!
! ===========================================================
!
!    POOL
!
! ===========================================================
!
      IF (KEYWIS ('POOL')) THEN
        CALL MSGERR(1,'Keyword POOL is not relevant anymore!')            40.31
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!     TEST       parameters for required test output
!
!     =============================================================
!
!                                  / -> IJ < [i] [j] > | < [k] >  \
!    TEST [itest] [itrace] POINTS <                                >  &
!                                  \    XY < [x] [y] >            /
!
!                          PAR 'fname'  S1D 'fname'  S2D 'fname'
!
!     ============================================================
!
      IF (KEYWIS ('TEST')) THEN
        CALL ININTG ('ITEST' , ITEST , 'STA', 30)
!       statements restructured:
        IF (ITEST.GE.30) THEN
          IF (ERRPTS.EQ.0.AND.IAMMASTER) THEN                             40.95 40.30 40.13
            ERRPTS = 16
            OPEN (ERRPTS, FILE='ERRPTS',
     &            STATUS='UNKNOWN', FORM='FORMATTED')
          ENDIF                                                           40.13
        ENDIF
        CALL ININTG ('ITRACE', ITRACE, 'UNC',  0)
        IF (ITRACE.GT.0) THEN
          LTRACE =.TRUE.
        ELSE
          LTRACE =.FALSE.
        ENDIF
        CALL INKEYW ('STA', ' ')
        IF (KEYWIS('POI')) THEN
          NPTST  = 0
          MPTST  = 50
          IPP    = 2                                                      40.80
          IF (OPTG.EQ.5) IPP = 1                                          40.80
          IF (.NOT.ALLOCATED(IARR)) ALLOCATE(IARR(IPP*MPTST))             40.80 40.31
!
          CALL RETSTP (IPP*MPTST, IARR, KGRPNT, KGRBND, XCGRID, YCGRID,   40.80 40.31
     &                 SPCSIG, SPCDIR)                                    40.31
          IF (STPNOW()) RETURN                                            34.01
          NPTSTA = MAX(1,NPTST)                                           40.00
          IF (.NOT.ALLOCATED(XYTST)) ALLOCATE(XYTST(IPP*NPTSTA))          40.80 40.31
          CALL SWCOPI (IARR,XYTST,IPP*NPTSTA)                             40.80 40.31
          DEALLOCATE(IARR)                                                40.31
        ENDIF
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     INTEst       parameters for required test output during input
!
!     ============================================================
!
!     INTE [intes]         (NOT documented)
!
!     ============================================================
!
      IF (KEYWIS ('INTE')) THEN
        CALL ININTG ('INTES' , INTES , 'STA', 30)
        GO TO 100
      ENDIF
!     ------------------------------------------------------------------
!     COTEst       parameters for required test output during computation
!
!     ============================================================
!
!     COTE [cotes]         (NOT documented)
!
!     ============================================================
      IF (KEYWIS ('COTE')) THEN
        CALL ININTG ('COTES' , ICOTES , 'STA', 30)
        GO TO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     OUTEst       parameters for required test output during output
!
!     ============================================================
!
!     OUTE [itest]         (NOT documented)
!
!     ============================================================
!
      IF (KEYWIS ('OUTE')) THEN
        CALL ININTG ('ITEST' , IOUTES , 'STA', 30)
        GO TO 100
      ENDIF

!     ============================================================

!      OUTPut OPTIons  'comment'  (TABle [field])  (BLOck  [ndec]  [len])   &

!      (SPEC  [ndec])

!     ============================================================

      IF (KEYWIS('OUTP')) THEN                                            40.13
        CALL IGNORE ('OPT')
        CALL INCSTR ('COMMENT', OUT_COMMENT, 'UNC', ' ')
        CALL INKEYW ('STA', ' ')
        IF (KEYWIS('TAB')) THEN
          CALL ININTG ('FIELD', FLD_TABLE, 'STA', 11)
          IF (FLD_TABLE.GT.16) CALL MSGERR (2, '[field] is too large')
          IF (FLD_TABLE.LT.8)  CALL MSGERR (2, '[field] is too small')
          WRITE (FLT_TABLE, 232) FLD_TABLE, FLD_TABLE-7
 232      FORMAT ('(E', I2, '.', I1, ')')
          IF (ITEST.GE.30) WRITE (PRINTF, 234) FLT_TABLE
 234      FORMAT (' Format floating point table: ', A)
          CALL INKEYW ('STA', ' ')
        ENDIF
        IF (KEYWIS('BLO')) THEN
          CALL ININTG ('NDEC', DEC_BLOCK, 'STA', 4)
          IF (DEC_BLOCK.GT.9) CALL MSGERR (2, '[ndec] is too large')
          CALL ININTG ('LEN', NLEN, 'STA', 200)
          IF (NLEN.GT.9999) CALL MSGERR (2, '[len] is too large')
          WRITE (FLT_BLOCK, 242) NLEN, DEC_BLOCK+7, DEC_BLOCK
 242      FORMAT ('(', I4, '(1X,E', I2, '.', I1, '))')
          IF (ITEST.GE.30) WRITE (PRINTF, 244) FLT_BLOCK
 244      FORMAT (' Format floating point block: ', A)
          CALL INKEYW ('STA', ' ')
        ENDIF
        IF (KEYWIS('SPEC')) THEN
          CALL ININTG ('NDEC', DEC_SPEC, 'STA', 5)
          IF (DEC_SPEC.GT.9) CALL MSGERR (2, '[ndec] is too large')
          CALL ININTG ('LEN', NLEN, 'STA', 200)
          IF (NLEN.GT.9999) CALL MSGERR (2, '[len] is too large')
          WRITE (FIX_SPEC, 252) NLEN, DEC_SPEC
 252      FORMAT ('(', I4, '(1X,I', I1, '))')
          IF (ITEST.GE.30) WRITE (PRINTF, 254) FIX_SPEC
 254      FORMAT (' Format spectral output: ', A)
          LENSPO = NLEN*(DEC_SPEC+1)
        ENDIF
        IF (KEYWIS('QUA')) THEN
          CALL MSGERR (3,
     &    'command OUTP QUANT is obsolete, use QUANTITY')
        ENDIF
        GOTO 100
      ENDIF                                                               40.13
!
!     ------------------------------------------------------------------
!
!     BOTTOM    definition of bottom grid
!
! ============================================================
!
! BOTtom ...         (OBSOLETE command)
!
! ============================================================
!
      IF (KEYWIS ('BOT')) THEN
        CALL MSGERR (2, 'command BOTTOM is replaced by INP BOT')
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     MODE  : Set STATionary, DYNamic (NONSTAtionary) or 1D SWAN model    32.02
!
! ========================================
!
!        | -> STAtionary    |     | -> TWODimensional |                   40.00
! MODE  <                    >   <                     >  (NOUPDATe)      40.07
!        |    NONSTationary |     |    ONEDimensional |                   40.00
!
! ========================================
!
      IF (KEYWIS ('MODE')) THEN                                           30.21
        CALL INKEYW ('STA',' ')
        IF (KEYWIS('NONST') .OR. KEYWIS ('DYN')) THEN                     30.70
          IF (NSTATM.EQ.0) CALL MSGERR (2, 'Mode Nonst incorrect here')   40.00
          NSTATM = 1                                                      40.00
          NSTATC = 1                                                      40.00
!
!         switch on flag for computation of default initial condition     30.75
          ICOND = 1                                                       30.75
!         no relaxation                                                   40.23
          PNUMS(30) = 0.                                                  40.23
        ELSEIF (KEYWIS ('STA')) THEN                                      40.00
          IF (NSTATM.EQ.1) CALL MSGERR (2, 'Mode STAT incorrect here')    40.00
          NSTATM = 0                                                      40.00
          CALL INKEYW ('STA',' ')                                         32.02
        ENDIF                                                             32.02
!
!       *** Logical ONED added for 1d-computations                        32.02
!
        CALL INKEYW ('STA', ' ')                                          40.00
        IF (KEYWIS ('ONED')) THEN                                         40.00
          ONED = .TRUE.                                                   32.02
          IF (PARLL) THEN                                                 40.30
             CALL MSGERR(4,'1D mode is not supported in parallel run')    40.30
             RETURN                                                       40.30
          END IF                                                          40.30
        ELSEIF (KEYWIS ('TWOD')) THEN
          ONED = .FALSE.                                                  32.02
        ENDIF
!
!       *** Logical ACUPDA added to avoid updating action densities       40.07
!
        CALL INKEYW ('STA', ' ')                                          40.07
        IF (KEYWIS ('NOUPDAT')) THEN                                      40.07
          ACUPDA = .FALSE.                                                40.07
        ENDIF
!
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
! =======================================================================
!
!            |    BSBT                      |
!   PROP    <                     | SEC |   |
!            |    GSE  [waveage] <  MIN  >  |
!            |                    |  HR |   |
!            |                    | DAY |   |
!
!            FLUXLIM
!
! =======================================================================
!
      IF (KEYWIS ('PROP')) THEN                                           40.02
        CALL INKEYW ('STA','    ')                                        40.02
        IF (KEYWIS ('BSBT')) THEN                                         40.02
          PROPSN = 1                                                      40.02
          PROPSS = 1                                                      40.02
        ELSE IF (KEYWIS ('GSE')) THEN                                     40.02
          IF (OPTG.NE.5 .AND. PROPSN.NE.3) THEN                           41.00 40.02
            CALL MSGERR(2,                                                40.02
     &      'Anti-GSE only allowed for S&L scheme.')                      40.02
          ENDIF                                                           40.02
          CALL ININTV('WAVEAGE', WAVAGE, 'STA', 0.)                       40.02
        ENDIF                                                             40.02
        IF (KEYWIS ('FLUXLIM')) THEN                                      40.23
           PROPFL   = 1                                                   40.23
           PNUMS(6) = 0.                                                  40.23
        END IF                                                            40.23
        GOTO 100                                                          40.02
      ENDIF                                                               40.02
!
!     ------------------------------------------------------------------
!
!     COORD : spherical or cartesian coordinates
!
! =============================================================
!
!   COORDinates  /  -> CARTesian               \   REPeating
!                \ SPHErical [rearth]   UM/QC  /
!
! =============================================================
!
      IF (KEYWIS ('COORD')) THEN                                          33.09
        CALL INKEYW ('STA',' ')                                           33.09
        IF (KEYWIS ('CART')) THEN                                         33.09
          KSPHER = 0                                                      33.09
        ELSE IF (KEYWIS ('SPHE')) THEN                                    33.09
          KSPHER = 1                                                      33.09
          CALL INREAL ('REARTH', REARTH, 'UNC', 0.)                       33.09
          LENDEG = REARTH * PI / 180.                                     33.09
!         change properties of output quantities Xp and Yp
          OVUNIT(1) = 'degr'
          OVLLIM(1) = -200.
          OVULIM(1) =  400.
          OVLEXP(1) = -180.
          OVHEXP(1) =  360.
          OVEXCV(1) = -999.
          OVUNIT(2) = 'degr'
          OVLLIM(2) = -100.
          OVULIM(2) =  100.
          OVLEXP(2) = -90.
          OVHEXP(2) =  90.
          OVEXCV(2) = -999.
          CALL INKEYW ('STA','CCM')                                       33.09
          IF (KEYWIS ('QC')) THEN                                         33.09
!           quasi-cartesian projection method
            PROJ_METHOD = 0                                               33.09
          ELSE IF (KEYWIS ('CCM')) THEN                                   33.09
!           uniform Mercator projection (default for spherical coordinates)
            PROJ_METHOD = 1                                               33.09
          ENDIF                                                           33.09
        ELSE                                                              33.09
          CALL WRNKEY
        ENDIF                                                             33.09
        CALL INKEYW ('STA',' ')                                           33.09
        IF (KEYWIS ('REP')) THEN                                          33.09
          KREPTX = 1                                                      33.09
        ENDIF                                                             33.09
        GO TO 100                                                         33.09
      ENDIF                                                               33.09
!
!     ------------------------------------------------------------------
!
!     **  Command COMPUTE   **                                            30.00
!
!       ==============================================================
!
!                  |  STATionary  [time]                      |
!       COMPute ( <                                            > )        40.00
!                  |                    | -> Sec  |           |
!                  |  ([tbegc] [deltc] <     MIn   > [tendc]) |
!                                       |    HR   |
!                                       |    DAy  |
!
!       ==============================================================
!
      IF (KEYWIS ('COMP')) THEN
        COMPUT = 'COMP'                                                   30.00
        RUNMADE = .TRUE.                                                  30.00
        CALL INKEYW ('STA','  ')
        IF (NSTATM.LE.0 .OR. KEYWIS('STAT')) THEN                         40.00
          IF (NSTATM.EQ.-1) NSTATM = 0                                    40.00
          IF (NSTATM.GT.0) CALL INCTIM (ITMOPT,'TIME',TINIC,'REQ',0.)     40.00
          IF (TINIC .LT. TIMCO) THEN                                      40.03
            print*,TINIC,TIMCO
            CALL MSGERR (2, '[time] before current time')                 40.03
            TINIC = TIMCO                                                 40.03
          ENDIF                                                           40.03
          TFINC = TINIC                                                   40.00
          TIMCO = TINIC                                                   40.00
          DT = 1.E10                                                      40.00
          RDTIM = 0.                                                      40.00
          NSTATC = 0                                                      40.00
          MTC = 1                                                         40.41
        ELSE
          CALL IGNORE ('NONST')
          IF (TIMCO .LT. -0.9E10) THEN                                    40.00
            CALL INCTIM (ITMOPT,'TBEGC',TINIC,'REQ',0.)                   40.03
          ELSE
            CALL INCTIM (ITMOPT,'TBEGC',TINIC,'STA',TIMCO)                40.03
          ENDIF
          IF (TINIC .LT. TIMCO) THEN                                      40.03
            CALL MSGERR (2, 'start time [tbegc] before current time')     40.03
            TINIC = TIMCO                                                 40.03
          ENDIF                                                           40.03
          CALL ININTV ('DELTC', DT, 'REQ', 0.)                            40.03
          CALL INCTIM (ITMOPT,'TENDC',TFINC,'REQ',0.)                     40.03
          NSTATC = 1                                                      40.00
!
!           *** tfinc must be greater than tinic **
          DIFF = TFINC - TINIC
          IF (DIFF .LE. 0.) CALL MSGERR (3,
     &        'start time [tbegc] greater or equal end time [tendc]')
!
!           **The number of computational steps is calculated
          RDTIM = 1./DT
          MTC = NINT ((TFINC - TINIC)/DT)                                 30.50
          IF (MOD(TFINC-TINIC,DT).GT.0.01*DT .AND.
     &            MOD(TFINC-TINIC,DT).LT.0.99*DT)
     &           CALL MSGERR (1,
     &           'DT is not a fraction of the computational period')
          TIMCO = TINIC
        ENDIF
        IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)                   40.00
        NCOMPT = NCOMPT + 1                                               40.41
        IF (NCOMPT.GT.50) CALL MSGERR (2,                                 40.41
     &                   'No more than 50 COMPUTE commands are allowed')  40.41
        RCOMPT(NCOMPT,1) = REAL(NSTATC)                                   40.41
        RCOMPT(NCOMPT,2) = REAL(MTC)                                      40.41
        RCOMPT(NCOMPT,3) = TFINC                                          40.41
        RCOMPT(NCOMPT,4) = TINIC                                          40.41
        RCOMPT(NCOMPT,5) = DT                                             40.41
!       set ITERMX equal to MXITST in case of stationary computations
!       and to MXITNS otherwise
        IF (NSTATC.EQ.0) THEN
          ITERMX = MXITST                                                 40.03
        ELSE
          ITERMX = MXITNS                                                 40.03
        ENDIF
        RETURN                                                            30.00
      ENDIF
!
!     ------------------------------------------------------------------
!
!     *** OBSTACLE   Definition of obstacles in comp grid. ***            30.61
!
! ============================================================
!
!             |  TRANSm [trcoef]                          |
! OBSTacle   <                                            |
!             |       | -> GODA [hgt] [alpha] [beta]       >  &
!             |  DAM <                                    |
!                     |    DANGremond [hgt] [slope] [Bk]  |
!                                                                         40.18
!                        | -> RSPEC        |                              40.18
!      ( REFLec [reflc] <                   > )   &
!                        |    RDIFF [pown] |                              40.13 40.18
!
!             LINe < [xp] [yp] >                                          40.18 40.09
!
! ============================================================
!
      IF (KEYWIS ('OBST')) THEN
!
        ALLOCATE(OBSTMP)                                                  40.31
        OBSTMP%TRCOEF(1) = 0.                                             40.31
        OBSTMP%TRCOEF(2) = 0.                                             40.31
        OBSTMP%TRCOEF(3) = 0.                                             40.31
        OBSTMP%RFCOEF(1) = 0.                                             40.31
        OBSTMP%RFCOEF(2) = 0.                                             40.31
        OBSTMP%RFCOEF(3) = 0.                                             40.31
        OBSTMP%RFCOEF(4) = 0.                                             40.31
        OBSTMP%RFCOEF(5) = 0.                                             40.31
        OBSTMP%RFCOEF(6) = 0.                                             40.31
!
!       data concerning transmission of energy over/through the obstacle
!
        CALL INKEYW ('STA', '  ')
        IF (KEYWIS ('TRANS')) THEN
          ITRAS = 0
          CALL INREAL ('TRCOEF', TRCF, 'REQ', 0.)
          IF ((TRCF.LT.0.) .OR. (TRCF.GT.1.)) THEN                        40.14
            CALL MSGERR(3,'Transmission coeff. [trans] is not allowed ')  40.14
            CALL MSGERR(3,'to be greater than 1 or smaller than 0!  ')    40.14
          ENDIF                                                           40.14
          OBSTMP%TRCOEF(1) = TRCF                                         40.31
        ELSE IF (KEYWIS ('DAM')) THEN
          CALL INKEYW ('STA', 'GODA')                                     40.66
          IF (KEYWIS ('DANG')) THEN                                       40.66
             ITRAS = 2                                                    40.66
             CALL INREAL ('HGT'  , HGT , 'REQ', 0.)                       40.66
             CALL INREAL ('SLOPE', SLP , 'REQ', 0.)                       40.66
             CALL INREAL ('BK'   , BK  , 'REQ', 0.)                       40.66
             OBSTMP%TRCOEF(1) = HGT                                       40.66
             OBSTMP%TRCOEF(2) = SLP                                       40.66
             OBSTMP%TRCOEF(3) = BK                                        40.66
          ELSE                                                            40.66
             CALL IGNORE ('GODA')                                         40.66
             ITRAS = 1
             CALL INREAL ('HGT'  , HGT , 'REQ', 0.  )                     30.80
             CALL INREAL ('ALPHA', OGAM, 'STA', 2.6 )                     30.60
             CALL INREAL ('BETA' , OBET, 'STA', 0.15)                     30.60
             OBSTMP%TRCOEF(1) = HGT                                       40.31
             OBSTMP%TRCOEF(2) = OGAM                                      40.31
             OBSTMP%TRCOEF(3) = OBET                                      40.31
          END IF                                                          40.66
        ELSE
!         if no transmission options are activated, there will be 0 transmission
          ITRAS = 0
          TRCF  = 0.                                                      40.41
          OBSTMP%TRCOEF(1) = TRCF                                         40.51
        ENDIF
        OBSTMP%TRTYPE = ITRAS                                             40.31
!                                                                         40.09
!       data on reflection by the obstacle, activated by keyword REFL     40.09
!
        CALL INKEYW ('REQ', '  ')                                         40.09
        IF (KEYWIS ('REFL')) THEN                                         40.09
          LREF = 1                                                        40.31
          IF (.NOT.FULCIR) THEN                                           40.09
            CALL MSGERR(3,'Reflections will only be calculated if     ')  40.09
            CALL MSGERR(3,'the spectral directions cover the full     ')  40.09
            CALL MSGERR(3,'circle.                                    ')  40.09
          ENDIF                                                           40.09
          CALL INREAL ('REFLC', REF0, 'STA', 1.)                          40.09
          DUM = REF0*REF0+TRCF*TRCF                                       40.51 40.14
          IF (ITRAS.EQ.0 .AND. DUM.GT.1.) THEN                            40.51 40.14
            CALL MSGERR(3,'Kt^2 + Kr^2 > 1 ')                             40.51 40.14
          ENDIF                                                           40.14
          OBSTMP%RFCOEF(1) = REF0                                         40.31
!                                                                         40.18
          CALL INKEYW ('REQ', '  ')                                       40.18
          IF (KEYWIS('RDIFF')) THEN                                       40.18
            LREFDIFF = 1                                                  40.31 40.18
            CALL INREAL ('POWN', POWN, 'REQ', 1.)                         40.18
            DUM = MOD(POWN,1.)                                            40.18
            IF (POWN.LT.0.) THEN                                          40.38
              CALL MSGERR(3,'Power POWN is not a positive number! ')      40.38
            ENDIF                                                         40.38
            IF (DUM.NE.0.) THEN                                           40.18
              CALL MSGERR(3,'Power POWN is not an integer number! ')      40.18
            ENDIF                                                         40.18
            OBSTMP%RFCOEF(2) = POWN                                       40.31
!                                                                         40.18
          ELSE                                                            40.18
            CALL IGNORE ('RSPEC')                                         40.18
            LREFDIFF = 0                                                  40.31 40.18
          ENDIF                                                           40.18
          OBSTMP%RFTYP2 = LREFDIFF                                        40.31
!
          CALL INKEYW ('REQ', '  ')                                       40.28
          CALL INKEYW ('STA', 'RFD')                                      40.28
          IF (KEYWIS('RFD')) THEN                                         40.28
            LRFRD = 1                                                     40.31 40.28
            CALL INREAL ('FD1', FD1, 'REQ', 0.565)                        40.28
            CALL INREAL ('FD2', FD2, 'REQ', 0.75)                         40.28
            CALL INREAL ('FD3', FD3, 'REQ', -21.)                         40.28
            CALL INREAL ('FD4', FD4, 'REQ', 0.066)                        40.28
            IF (FD3.GE.0.) THEN                                           40.28
              CALL MSGERR(3,'Positive frequency dependency! ')            40.28
            ENDIF                                                         40.28
            OBSTMP%RFCOEF(3) = FD1                                        40.31
            OBSTMP%RFCOEF(4) = FD2                                        40.31
            OBSTMP%RFCOEF(5) = FD3                                        40.31
            OBSTMP%RFCOEF(6) = FD4                                        40.31
          ELSE                                                            40.28
            LRFRD = 0                                                     40.31 40.28
          ENDIF                                                           40.28
          OBSTMP%RFTYP3 = LRFRD                                           40.31
!                                                                         40.18
          CALL INKEYW ('REQ', '  ')                                       40.09
!                                                                         40.09
          IF ((REF0.LT.0.) .OR. (REF0.GT.1.)) THEN                        40.09
            CALL MSGERR(3,'Reflection coeff. [reflc] is not allowed ')    40.09
            CALL MSGERR(3,'to be greater than 1 or smaller than 0!  ')    40.09
          ENDIF                                                           40.09
        ELSE
!         if there is no keyword REFL, there will be no reflection        40.09
          LREF = 0
        ENDIF                                                             40.09
        OBSTMP%RFTYP1 = LREF                                              40.31
!
!       location of obstacles
!
!       *** NUMCOR : Number of corners ***
        NUMCOR = 0
        IF (KEYWIS ('LIN')) THEN
          FRST%X = 0.                                                     40.31
          FRST%Y = 0.                                                     40.31
          NULLIFY(FRST%NEXTXY)                                            40.31
          CURR => FRST                                                    40.31
          DO
!           read coordinates of one corner point of the obstacle
            CALL READXY ('XP', 'YP', XP, YP, 'REP', -1.E10, -1.E10)
            IF (XP.LT.-.9E10) GOTO 101
            NUMCOR = NUMCOR + 1
            ALLOCATE(TMP)                                                 40.31
            TMP%X = XP                                                    40.31
            TMP%Y = YP                                                    40.31
            NULLIFY(TMP%NEXTXY)                                           40.31
            CURR%NEXTXY => TMP                                            40.31
            CURR => TMP                                                   40.31
          ENDDO
        ENDIF
 101    OBSTMP%NCRPTS = NUMCOR                                            40.31
!       store coordinates for corner points in array of obstacle data     40.31
        ALLOCATE(OBSTMP%XCRP(NUMCOR))                                     40.31
        ALLOCATE(OBSTMP%YCRP(NUMCOR))                                     40.31
        CURR => FRST%NEXTXY                                               40.31
        DO JJ = 1, NUMCOR                                                 40.31
            OBSTMP%XCRP(JJ) = CURR%X                                      40.31
            OBSTMP%YCRP(JJ) = CURR%Y                                      40.31
            CURR => CURR%NEXTXY                                           40.31
        END DO                                                            40.31
        DEALLOCATE(TMP)                                                   40.31
        NULLIFY(OBSTMP%NEXTOBST)                                          40.31
        IF (NUMCOR .LE. 1) THEN
          CALL MSGERR(1,'No corner points for obstacle were found')
        ELSE
!         *** NUMOBS : Number of obstacles ***
          NUMOBS = NUMOBS + 1
          IF ( .NOT.LOBST ) THEN                                          40.31
             FOBSTAC = OBSTMP                                             40.31
             COBST => FOBSTAC                                             40.31
             LOBST = .TRUE.                                               40.31
          ELSE                                                            40.31
             COBST%NEXTOBST => OBSTMP                                     40.31
             COBST => OBSTMP                                              40.31
          END IF                                                          40.31
        ENDIF
!
        GO TO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     ***'INITial conditions'  Definition of initial conditions  ***
!     *** for MODE DYNAMIC
!
      IF (KEYWIS ('INIT')) THEN                                           30.61
        CALL INITVA( AC2, SPCSIG, SPCDIR, KGRPNT )                        40.80 40.31
        IF (STPNOW()) RETURN                                              34.01
        GO TO 100
      ENDIF
!
!     ------------------------------------------------------------------
!     HOTFile    write current wave field to file for future use as
!     initial cond.
      IF (KEYWIS('REST') .OR. KEYWIS ('BACK') .OR. KEYWIS('HOTF')         40.00
     &    .OR. KEYWIS('SAVE')) THEN                                       40.00
        IF (MXC.LE.0 .AND. OPTG.NE.5) THEN                                40.80
          CALL MSGERR (2, 'command CGRID must precede this command')
        ELSEIF (MCGRD .LE. 1 .AND. nverts .LE. 0) THEN                    40.80
          CALL MSGERR(2,                                                  40.80
     &    ' command READ BOT or READ UNSTRUC must precede this command')  40.80
        ELSE
          CALL BACKUP( AC2, SPCSIG, SPCDIR, KGRPNT, XCGRID, YCGRID )      40.31 30.90
          IF (STPNOW()) RETURN                                            34.01
        ENDIF
        GO TO 100
      ENDIF
!     ------------------------------------------------------------------
!
!     INPUT    definition of input grids
!
! ============================================================================
!
!   INPgrid                                                                &
!      BOTtom / WLEVel / CURrent / VX / VY / FRiction / WInd / WX / WY     &
!      NPLAnts                                                             &
!      | REG [xpinp] [ypinp] [alpinp]  [mxinp] [myinp]  [dxinp] [dyinp] |
!     <  CURVilinear [stagrx] [stagry] [mxinp] [myinp]                   > &
!      | UNSTRUCtured                                                   |
!      (NONSTATionary [tbeginp] [deltinp] SEC/MIN/HR/DAY [tendinp])
!
! ============================================================================
!
      IF (KEYWIS ('INP')) THEN
        CALL IGNORE ('GRID')
        CALL INKEYW ('STA', ' ')
        IGR2 = 0                                                          10.26
        IF (KEYWIS ('BOT')) THEN
          IGRD = 1
          PSNAME = 'BOTTGRID'
        ELSE IF (KEYWIS ('CUR')) THEN
          IGRD = 2
          IGR2 = 3                                                        10.26
          PSNAME = 'VXGRID  '
        ELSE IF (KEYWIS ('VX')) THEN
          IGRD = 2
          PSNAME = 'VXGRID  '
        ELSE IF (KEYWIS ('VY')) THEN
          IGRD = 3
          PSNAME = 'VYGRID  '
        ELSE IF (KEYWIS ('FR')) THEN
          IGRD = 4
          PSNAME = 'FRICGRID'
        ELSE IF (KEYWIS ('WI')) THEN
          IGRD = 5
          IGR2 = 6                                                        10.26
          PSNAME = 'WXGRID  '
        ELSE IF (KEYWIS ('WX')) THEN
          IGRD = 5
          PSNAME = 'WXGRID  '
        ELSE IF (KEYWIS ('WY')) THEN
          IGRD = 6
          PSNAME = 'WYGRID  '
        ELSE IF (KEYWIS ('WLEV')) THEN                                    20.38
          IGRD = 7
          PSNAME = 'WLEVGRID'
!       note: 8 and 9 are for coordinates                                 40.03
        ELSE IF (KEYWIS ('ASTD')) THEN                                    40.03
!         air-sea temperature difference                                  40.03
          IGRD = 10
          PSNAME = 'ASTDGRID'
        ELSE IF (KEYWIS ('NPLA')) THEN                                    40.55
!         number of plants per square meter                               40.55
          IGRD = 11                                                       40.55
          PSNAME = 'NPLAGRID'                                             40.55
        ELSE
          IGRD = 1
          PSNAME = 'BOTTGRID'
        ENDIF
!
        CALL SINPGR (IGRD, IGR2, PSNAME)                                  40.31
        IF (STPNOW()) RETURN                                              34.01
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!     READ   reading depths, coordinates and/or currents
!
!   ============================================================================
!
!   READinp    BOTtom/WLevel/CURrent/FRiction/WInd/COORdinates               &
!              NPLAnts                                                       &
!        [fac]  / 'fname1'        \
!               \ SERIES 'fname2' /  [idla] [nhedf] ([nhedt]) (nhedvec])     &
!        FREE / FORMAT 'form' / [idfm] / UNFORMATTED
!
!   or
!
!                          | -> ADCirc
!   READgrid UNSTRUCtured <  TRIAngle \
!                          | EASYmesh / 'fname'
!
!   ============================================================================
!
      IF (KEYWIS ('READ')) THEN
!
!       --- check whether computational grid has been defined             40.80
!
        IF ( .NOT.LOGCOM(2) ) THEN                                        40.80
           CALL MSGERR (3, '* define computational grid before      *')   40.80
           CALL MSGERR (3, '* reading coordinates or input grids    *')   40.80
        ENDIF                                                             40.80
!
        CALL INKEYW ('REQ',' ')
        IF (KEYWIS ('UNSTRUC')) THEN                                      40.80
!
           IF (OPTG.NE.5) THEN                                            40.80
              CALL MSGERR (3,                                             40.80
     &     'define computational grid before reading unstructured grid')  40.80
           ENDIF                                                          40.80
           IF (ONED) THEN                                                 40.80
              CALL MSGERR (4,                                             40.80
     &            '1D-simulation cannot be done with unstructured grid')  40.80
              RETURN                                                      40.80
           ENDIF                                                          40.80
           IF (PARLL) THEN                                                40.80
              CALL MSGERR(4,
     &             'Unstructured grid is not supported in parallel run')  40.80
              RETURN                                                      40.80
           ENDIF                                                          40.80
!
!          --- read unstructured grid                                     40.80
!
           LOGCOM(5) = .TRUE.                                             40.80
           CALL INKEYW ('STA', 'ADC')                                     40.80
           IF (KEYWIS('ADC')) THEN                                        40.80
              grid_generator = meth_adcirc                                40.80
!             bottom topography will be taken from fort.14, so            40.80
              LOGCOM(3) = .TRUE.                                          40.80
              IGTYPE(1) = 3                                               40.80
              LEDS(1)   = 2                                               40.80
           ELSEIF (KEYWIS('TRIA')) THEN                                   40.80
              grid_generator = meth_triangle                              40.80
              CALL INCSTR ('FNAME',FILENM,'REQ', ' ')                     40.80
           ELSEIF (KEYWIS('EASY')) THEN                                   40.80
              grid_generator = meth_easy                                  40.80
              CALL INCSTR ('FNAME',FILENM,'REQ', ' ')                     40.80
           ENDIF                                                          40.80
           CALL SwanReadGrid ( FILENM, LENFNM )                           40.80
           IF (STPNOW()) RETURN                                           40.80
           GOTO 444                                                       40.80
        ENDIF
!
        CALL SREDEP ( LWINDR, LWINDM, LOGCOM )                            40.31 40.02
        IF (STPNOW()) RETURN                                              34.01
!       --- call CGINIT once as soon as both coordinates and bottom       40.31
!           have been read and AC2 is not allocated                       40.31
 444    IF (OPTG .EQ. 1) THEN                                             40.80
!         regular grid
          IF (LOGCOM(3) .AND. .NOT. LOGCOM(6)) THEN                       30.90
            CALL CGINIT(LOGCOM)                                           40.31 30.90
            IF (STPNOW()) RETURN                                          34.01
          ENDIF                                                           30.90
        ELSEIF (OPTG.EQ.3) THEN                                           40.80
!         cuvilinear grid
          IF (LOGCOM(3) .AND. .NOT. LOGCOM(4)) THEN
            CALL MSGERR (3, '** Give CGRID command and               *')  40.00
            CALL MSGERR (3, '** read curvilinear coordinates         *')
            CALL MSGERR (3, '** before reading the bottom grid       *')
          ELSE IF (LOGCOM(2) .AND. LOGCOM(3) .AND.
     &             LOGCOM(4) .AND. .NOT. LOGCOM(6)) THEN
            CALL CGINIT(LOGCOM)                                           40.31 30.90
            IF (STPNOW()) RETURN                                          34.01
          ENDIF
        ELSEIF (OPTG.EQ.5) THEN                                           40.80
!         unstructured grid
!
          IF ( .NOT.LOGCOM(5) ) THEN                                      40.80
            CALL MSGERR (3, '* read unstructured grid by means of    *')  40.80
            CALL MSGERR (3, '* SMS/ADCIRC, Triangle or Easymesh      *')  40.80
            IF ( LOGCOM(4) )                                              40.80
     &      CALL MSGERR (3, '* instead of curvilinear coordinates    *')  40.80
          ELSEIF ( LOGCOM(5) .AND. LOGCOM(2) .AND.                        40.80
     &            .NOT.LOGCOM(4) .AND. .NOT.LOGCOM(6) ) THEN              40.80
            CALL SwanInitCompGrid (LOGCOM)                                40.80
            IF (STPNOW()) RETURN                                          40.80
            CALL SwanCreateEdges                                          40.80
            IF (STPNOW()) RETURN                                          40.80
            CALL SwanGridTopology                                         40.80
            IF (STPNOW()) RETURN                                          40.80
!
!           --- the computational grid is included in output data         40.80
            IF ( nverts.GT.0 ) THEN                                       40.80
               ALLOCATE(OPSTMP)                                           40.80
               OPSTMP%PSNAME = 'COMPGRID'                                 40.80
               OPSTMP%PSTYPE = 'U'                                        40.80
               OPSTMP%MIP    = nverts                                     40.80
               ALLOCATE(OPSTMP%XP(nverts))                                40.80
               ALLOCATE(OPSTMP%YP(nverts))                                40.80
               DO JJ = 1, nverts                                          40.80
                  OPSTMP%XP(JJ) = xcugrd(JJ)                              40.80
                  OPSTMP%YP(JJ) = ycugrd(JJ)                              40.80
               ENDDO                                                      40.80
               NULLIFY(OPSTMP%NEXTOPS)                                    40.80
               IF ( .NOT.LOPS ) THEN                                      40.80
                  FOPS = OPSTMP                                           40.80
                  COPS => FOPS                                            40.80
                  LOPS = .TRUE.                                           40.80
               ELSE                                                       40.80
                  COPS%NEXTOPS => OPSTMP                                  40.80
                  COPS => OPSTMP                                          40.80
               ENDIF                                                      40.80
            ENDIF                                                         40.80
          ENDIF                                                           40.80
        ENDIF
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     CGRID     definition of computational grid
!
! ===========================================================================
!
!          / REGular [xpc] [ypc] [alpc] [xlenc] [ylenc] [mxc] [myc] \
!   CGRID <  CURVilinear [mxc] [myc]  [excval]                       >   &
!          \ UNSTRUCtured                                           /
!
!          / CIRcle               \
!          \ SECtor [dir1] [dir2] /  [mdc]  [flow]  [fhig]  [msc]
!
! ===========================================================================
!
      IF (KEYWIS('CGRID') .OR. KEYWIS('GRID')) THEN                       30.20
!       ver 20.67: command reorganised in view of cycle 2 (parametric &
!                  curvilinear)
!       ver 30.20: order of data changed
!       [xpc], [ypc], [alpc] moved                                        30.20
!
!
!       ver 30.21: CURVILINEAR                                            30.21
!
!
        LOGCOM(2) = .TRUE.                                                30.21
        CALL INKEYW ('STA', 'REG')                                        30.60
        IF (KEYWIS('CURV')) THEN                                          30.21
          OPTG  = 3                                                       30.21
          XPC   = 0.
          YPC   = 0.
          ALPC  = 0.                                                      30.21
          XCLEN = 0.
          YCLEN = 0.
          IF (ITEST.GE.30) THEN                                           40.41
          CALL MSGERR(1,'there is an unresolved problem with the')        40.08
          CALL MSGERR(1,'curvilinear mode, so use with caution.')         40.08
          CALL MSGERR(1,'This problem occurs when the upwind point')      40.08
          CALL MSGERR(1,'has not been solved for yet, because it')        40.08
          CALL MSGERR(1,'falls into a different sweep. Normally,')        40.08
          CALL MSGERR(1,'an upwind point will always fall within')        40.08
          CALL MSGERR(1,'the same sweep, but this is not necessarily')    40.08
          CALL MSGERR(1,'the case when using curvilinear coordinates.')   40.08
          CALL MSGERR(1,'To state the problem another way:')              40.08
          CALL MSGERR(1,'the axis of the grid bends about the')           40.08
          CALL MSGERR(1,'direction associated with a particular')         40.08
          CALL MSGERR(1,'directional bin.')                               40.08
          END IF                                                          40.41
        ELSEIF (KEYWIS('UNSTRUC')) THEN                                   40.80
          OPTG  = 5                                                       40.80
          XPC   = 0.                                                      40.80
          YPC   = 0.                                                      40.80
          ALPC  = 0.                                                      40.80
          XCLEN = 0.                                                      40.80
          YCLEN = 0.                                                      40.80
        ELSE
          CALL IGNORE ('REG')                                             20.67
          OPTG = 1
          IF (KSPHER.EQ.0) THEN                                           40.13
            CALL READXY ('XPC', 'YPC', XPC, YPC, 'STA', 0., 0.)
            CALL INREAL ('ALPC',ALPC,'UNC',0.)
          ELSE                                                            40.13
!           spherical coordinates; [xpc] and [ypc] are required           40.13
            CALL READXY ('XPC', 'YPC', XPC, YPC, 'REQ', 0., 0.)           40.13
            CALL INREAL ('ALPC',ALPC,'STA',0.)                            40.13
          ENDIF                                                           40.13
          CALL INREAL('XLENC',XCLEN,'RQI',0.)                             30.20
          IF (ONED) THEN                                                  32.02
            CALL INREAL('YLENC',YCLEN,'STA',0.)                           40.13
            IF (YCLEN .NE. 0) THEN                                        32.02
              CALL MSGERR (1, '1D-simulation: [ylenc] set to zero !')     32.02
            ENDIF                                                         32.02
            YCLEN = 0.                                                    40.00
          ELSE
            CALL INREAL('YLENC',YCLEN,'RQI',0.)                           30.20
          ENDIF
!
!         ALPC is made to be between -PI and PI
!
          ALTMP = ALPC / 360.
          ALPC = PI2 * (ALTMP - NINT(ALTMP))
          CVLEFT = .TRUE.                                                 40.00
        ENDIF                                                             30.21
!PUN!
!PUN        IF (OPTG.NE.5) THEN                                               40.95
!PUN           CALL MSGERR(4,
!PUN     &               'Structured grid is not supported in parallel run')  40.95
!PUN           RETURN                                                         40.95
!PUN        ENDIF                                                             40.95
!
        IF (OPTG.EQ.5) GOTO 555                                           40.80
!
!       ***** MXC is the number of steps in X direction *****
!       ***** MYC is the number of steps in Y direction *****
!
        CALL ININTG('MXC',MXS,'RQI',0)
!
        IF (ONED) THEN                                                    32.02
          CALL ININTG('MYC',MYS,'STA',0)                                  40.13
          IF (MYS .NE. 0) THEN                                            32.02
            CALL MSGERR (1, '1D-simulation: [myc] set to zero !')         32.02
          ENDIF                                                           32.02
          MYS = 0                                                         32.02
        ELSE                                                              40.13
          CALL ININTG('MYC',MYS,'RQI',-1)
        ENDIF                                                             32.02
!
        IF (KREPTX.EQ.1) THEN                                             33.09
          MXC = MXS                                                       33.09
        ELSE
          MXC = MXS+1
        ENDIF
        MYC = MYS+1
        MMCGR = MXC*MYC
        DX  = XCLEN/MXS
!                                                                         32.02
        IF (ONED) THEN                                                    32.02
          DY  = DX                                                        32.02
        ELSE                                                              32.02
          DY  = YCLEN/MYS
        ENDIF                                                             32.02
 555    CONTINUE                                                          40.80
!
!       for curvilinear grid, read exception values for grid point coordinates
!
        IF (OPTG.EQ.3) THEN                                               30.60
          CALL INKEYW ('STA', ' ')                                        30.60
          IF (KEYWIS ('EXC')) THEN                                        30.60
            CALL INREAL ('EXCVAL', EXCFLD(8), 'REQ', 0.)                  30.60
            CALL INREAL ('EXCVAL', EXCFLD(9), 'STA', EXCFLD(8))           30.60
          ENDIF                                                           30.60
        ENDIF                                                             30.60
!
        CALL INKEYW ('STA', 'CIR')                                        20.67
        IF (KEYWIS('SEC')) THEN
          FULCIR = .FALSE.
          CALL INREAL ('DIR1', SPDIR1, 'REQ', 0.)                         20.43
          CALL INREAL ('DIR2', SPDIR2, 'REQ', 0.)                         20.56
          CALL ININTG('MDC',MDC,'RQI',0)
        ELSE                                                              20.43
          CALL IGNORE ('CIR')                                             20.67
          FULCIR = .TRUE.
          CALL ININTG('MDC',MDC,'RQI',0)
        ENDIF
!
!       ***** MSC is the number of frequencies in sigma direction (logaritmic)
!
        CALL INREAL('FLOW',FRLOW,'STA',-999.)                             40.31
        CALL INREAL('FHIGH',FRHIG,'STA',-999.)                            40.31
        CALL ININTG('MSC',MSS,'STA',-999)                                 40.31
        IVAL = MAX(ABS(MSS*NINT(FRLOW)),ABS(MSS*NINT(FRHIG)),             40.31
     &             ABS(NINT(FRLOW)*NINT(FRHIG)))                          40.31
        IF (IVAL.GE.999**2) THEN                                          40.31
           CALL MSGERR(3,                                                 40.31
     &                 'At least, FLOW and FHIGH or MSC must be given!')  40.31
        END IF                                                            40.31
        IVAL = MAX(ABS(NINT(FRLOW)),ABS(NINT(FRHIG)),ABS(MSS))            40.31
        IF (IVAL.NE.999) THEN                                             40.31
           GAMMA = EXP(ALOG(FRHIG/FRLOW)/REAL(MSS)) - 1.                  40.41
           WRITE (PRINTF,'(A,(F7.4))')                                    40.41
     &                       ' Resolution in sigma-space: df/f = ',GAMMA  40.41
        ELSE IF (MSS.EQ.-999) THEN                                        40.31
           MSS = NINT(ALOG(FRHIG/FRLOW)/ALOG(1.1))                        40.31
           WRITE (PRINTF,'(A,I3)')                                        40.31
     &                  ' Number of meshes in sigma-space: MSC-1 = ',MSS  40.41 40.31
        ELSE IF (FRLOW.EQ.-999.) THEN                                     40.31
           FRLOW = FRHIG/EXP(ALOG(1.1)*REAL(MSS))                         40.31
           WRITE (PRINTF,'(A,(F7.4))')                                    40.31
     &               ' Lowest discrete frequency (in Hz): FLOW = ',FRLOW  40.31
        ELSE IF (FRHIG.EQ.-999.) THEN                                     40.31
           FRHIG = FRLOW*EXP(ALOG(1.1)*REAL(MSS))                         40.31
           WRITE (PRINTF,'(A,(F7.4))')                                    40.31
     &             ' Highest discrete frequency (in Hz): FHIGH = ',FRHIG  40.31
        END IF                                                            40.31
        SLOW = 2.*PI*FRLOW
        SHIG = 2.*PI*FRHIG
        MSC  = MSS+1
!
!       ***** MDC is the number of steps in theta direction as part of a circle
!
        IF (FULCIR) THEN
          DDIR  = PI2 / MDC                                               20.43
!
!         modification of SPDIR1 first installed with version 30.72, then
!         reversed and reinstalled with 40.13
!         purpose: prevent problem with grids under 45 degrees            40.13
          SPDIR1 = ALPC + 0.5 * DDIR                                      30.72 40.13
        ELSE
          IF (BNAUT) THEN                                                 30.70
!           swap values of SPDIR1 and SPDIR2, and transform               30.70
            TMPDIR = SPDIR1                                               30.70
            SPDIR1 = 180. + DNORTH - SPDIR2                               30.70
            SPDIR2 = 180. + DNORTH - TMPDIR                               30.70
          ENDIF                                                           30.70
          SPDIR1 = SPDIR1 * PI / 180.                                     30.50
          SPDIR2 = SPDIR2 * PI / 180.                                     30.50
          IF (SPDIR2.LT.SPDIR1) SPDIR2 = SPDIR2 + PI2                     20.56
          DDIR = (SPDIR2-SPDIR1) / REAL(MDC)                              20.56
          MDC = MDC + 1                                                   20.56
        ENDIF
!
        IF(.NOT.ALLOCATED(SPCSIG)) ALLOCATE(SPCSIG(MSC))                  40.31
        IF(.NOT.ALLOCATED(SPCDIR)) ALLOCATE(SPCDIR(MDC,6))                40.31
        CALL SSFILL(SPCSIG,SPCDIR)                                        40.31 30.90
!
        IF (ITEST.GE. 20) THEN
          IF(OPTG .EQ. 1)WRITE (PRINTF,6048)                              30.21
          IF(OPTG .EQ. 3)WRITE (PRINTF,6049)                              30.21
          IF(OPTG .EQ. 5)WRITE (PRINTF,6050)                              40.80
          WRITE (PRINTF,6045) SLOW, SHIG, FRINTF
          IF (OPTG.NE.5) WRITE (PRINTF,6046) MXC,MYC,MDC,MSC
          IF (OPTG.NE.5) WRITE (PRINTF,6047) DX,DY,DDIR
          IF (OPTG.EQ.5) WRITE (PRINTF,6051) MDC,MSC,DDIR                 40.80
 6048     FORMAT ('GRID: REGULAR RECTANGULAR')                            30.21
 6049     FORMAT ('GRID: CURVILINEAR')                                    30.21
 6050     FORMAT ('GRID: UNSTRUCTURED')                                   40.80
 6045     FORMAT (' S-low: ', F6.3,' S-hig: ', F6.3, ' frintf: ', F6.3)
 6046     FORMAT (' MXC: ',I6,' MYC: ',I6,' MDC: ',I6,' MSC: ',I6)
 6047     FORMAT (' DX: ',E12.4,' DY: ',E12.4, ' DDIR: ', F6.3)
 6051     FORMAT (' MDC: ',I6,' MSC: ',I6, ' DDIR: ', F6.3)               40.80
        ENDIF
!
        IF(.NOT.ALLOCATED(XCGRID)) ALLOCATE(XCGRID(MXC,MYC))              40.31
        IF(.NOT.ALLOCATED(YCGRID)) ALLOCATE(YCGRID(MXC,MYC))              40.31
!
        IF (OPTG .EQ. 1) THEN                                             30.60
!         *** The coordinates of computational points in ***
!         *** regular grid are computed                  ***
          COSPC = COS(ALPC)                                               7/JAN
          SINPC = SIN(ALPC)                                               7/JAN
          DO 200 J = 1,MYC
            DO 205 I = 1, MXC
              VALX = XPC + COSPC*(I-1)*DX - SINPC*(J-1)*DY                7/JAN
              VALY = YPC + SINPC*(I-1)*DX + COSPC*(J-1)*DY                7/JAN
              XCGRID(I,J)=VALX                                            40.31
              YCGRID(I,J)=VALY                                            40.31
 205        CONTINUE
 200      CONTINUE
        ENDIF
!
!       *** the computational grid is included in output data  ***        40.31
        IF ((MXC.GT.0) .AND. (MYC.GT.0) ) THEN                            40.31
           ALLOCATE(OPSTMP)                                               40.31
           OPSTMP%PSNAME = 'COMPGRID'                                     40.31
           IF (OPTG .EQ. 1) THEN                                          40.31
              OPSTMP%PSTYPE = 'F'                                         40.31
              OPSTMP%OPR(1) = XPC                                         40.31
              OPSTMP%OPR(2) = YPC                                         40.31
              OPSTMP%OPR(3) = XCLEN                                       40.31
              OPSTMP%OPR(4) = YCLEN                                       40.31
              OPSTMP%OPR(5) = ALPC                                        40.31
           ELSE IF (OPTG.EQ.3) THEN                                       40.80 40.31
              OPSTMP%PSTYPE = 'H'                                         40.31
              OPSTMP%OPR(1) = FLOAT(MXC-1)                                40.31
              OPSTMP%OPR(2) = FLOAT(MYC-1)                                40.31
              OPSTMP%OPR(3) = 0.                                          40.31
              OPSTMP%OPR(4) = 0.                                          40.31
              OPSTMP%OPR(5) = ALPC                                        40.31
           ENDIF                                                          40.31
           OPSTMP%OPI(1) = MXC                                            40.31
           OPSTMP%OPI(2) = MYC                                            40.31
           ALLOCATE(OPSTMP%XP(0))                                         40.31
           ALLOCATE(OPSTMP%YP(0))                                         40.31
           NULLIFY(OPSTMP%NEXTOPS)                                        40.31
           IF ( .NOT.LOPS ) THEN                                          40.31
              FOPS = OPSTMP                                               40.31
              COPS => FOPS                                                40.31
              LOPS = .TRUE.                                               40.31
           ELSE                                                           40.31
              COPS%NEXTOPS => OPSTMP                                      40.31
              COPS => OPSTMP                                              40.31
           END IF                                                         40.31
        ENDIF
!
        GOTO 100
      ENDIF
!
!     -----------------------------------------------------------------------
!
!     BOUNdary  defining boundary conditions                              20.63
!
      IF (KEYWIS ('BOU')) THEN
        ITMP1  = MXC                                                      40.31
        ITMP2  = MYC                                                      40.31
        ITMP3  = MCGRD                                                    40.31
        ITMP4  = NGRBND                                                   40.31
        MXC    = MXCGL                                                    40.31
        MYC    = MYCGL                                                    40.31
        MCGRD  = MCGRDGL                                                  40.31
        NGRBND = NGRBGL                                                   40.31
        CALL SWBOUN ( XGRDGL, YGRDGL, KGRPGL, XYTST, KGRBGL )             40.31
        IF (STPNOW()) RETURN                                              34.01
        MXC    = ITMP1                                                    40.31
        MYC    = ITMP2                                                    40.31
        MCGRD  = ITMP3                                                    40.31
        NGRBND = ITMP4                                                    40.31
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     LIM     setting parameters in conjunction with action limiter
!
! ============================================================
!
!   LIMiter [ursell] [qb]                                                 40.16
!
! ============================================================
!
      IF (KEYWIS ('LIM')) THEN                                            40.16
        CALL INREAL('URSELL',PTRIAD(3),'UNC',0.)                          40.16
        CALL INREAL('QB'    ,PNUMS(28),'UNC',0.)                          40.16
        GOTO 100                                                          40.16
      ENDIF                                                               40.16
!
!     ------------------------------------------------------------------
!
!     NUM       setting numerical parameters for schemes, solvers, etc.
!
! =====================================================================
!
!             | -> ACCUR [drel] [dhoval] [dtoval] [npnts]            |    40.41 40.03
!   NUMeric (<                                                        > & 40.41
!             | STOPC [dabs] [drel] [curvat] [npnts] [dtabs] [curvt] |    40.93 40.41
!
!                    | -> STAT  [mxitst] [alfa] |                         40.23
!                   <                            >  [limiter]   )     &   40.03
!                    | NONSTat  [mxitns]        |
!
!           ( DIRimpl [cdd] [cdlim]  WNUMber                       )  &
!
!           ( REFRLim [frlim] [power]          (NOT documented)    )  &
!
!           | -> SIGIMpl [css] [eps2] [outp] [niter]               |
!          (<                                                      >) &
!           |    SIGEXpl [css] [cfl]           (NOT documented)    |
!           |                                                      |
!           |    FIL     [diffc]               (NOT documented)    |
!
!
!           ( SETUP [eps2] [outp] [niter]                          )      30.82
!
! =====================================================================
!
      IF (KEYWIS ('NUM')) THEN
        CALL INKEYW ('REQ','  ')
!       *** accuracy and criterion to terminate the iteration ***
        IF (KEYWIS ('ACCUR')) THEN
          PNUMS(21) = 0.                                                  40.41
          CALL INREAL ('DREL'   , PNUMS(1) , 'UNC', 0.)
          CALL INREAL ('DHOVAL' , PNUMS(15), 'UNC', 0.)                   30.82
          CALL INREAL ('DTOVAL' , PNUMS(16), 'UNC', 0.)                   30.82
          CALL INREAL ('NPNTS'  , PNUMS(4) , 'UNC', 0.)
          CALL INKEYW ('STA', 'STAT')                                     40.03
          IF (KEYWIS ('STAT')) THEN                                       40.03
            CALL ININTG ('MXITST' , MXITST   , 'UNC', 0 )                 40.80 40.03
            CALL INREAL ('ALFA'   , PNUMS(30), 'UNC', 0.)                 40.23
          ELSE IF (KEYWIS ('ITERMX')) THEN                                40.03
            CALL ININTG ('MXITST' , MXITST   , 'REQ', 0 )                 40.03
            MXITNS = MXITST                                               40.03
            CALL MSGERR (1, '[itermx] is replaced; see user manual')      40.03
          ELSE IF (KEYWIS ('NONST')) THEN                                 40.03
            CALL ININTG ('MXITNS' , MXITNS   , 'REQ', 0 )                 40.03
          ENDIF
          CALL INREAL ('LIMITER', PNUMS(20), 'UNC', 0.)                   30.70
        ELSE IF (KEYWIS ('STOPC')) THEN                                   40.41
          PNUMS(21) = 1.                                                  40.41
          CALL INREAL ('DABS'   , PNUMS(2) , 'STA', 0.00)                 40.41
          CALL INREAL ('DREL'   , PNUMS(1) , 'STA', 0.01)                 40.41
          CALL INREAL ('CURVAT' , PNUMS(15), 'STA', 0.005)                40.41
          CALL INREAL ('NPNTS'  , PNUMS(4) , 'UNC', 0.)                   40.41
          CALL INREAL ('DTABS'  , PNUMS(3) , 'STA', 1000.)                40.93
          CALL INREAL ('CURVT'  , PNUMS(16), 'STA', 1000.)                40.93
          CALL INKEYW ('STA', 'STAT')                                     40.41
          IF (KEYWIS ('STAT')) THEN                                       40.41
            CALL ININTG ('MXITST' , MXITST   , 'STA', 50)                 40.41
            CALL INREAL ('ALFA'   , PNUMS(30), 'UNC', 0.)                 40.41
          ELSE IF (KEYWIS ('ITERMX')) THEN                                40.41
            CALL ININTG ('MXITST' , MXITST   , 'REQ', 0 )                 40.41
            MXITNS = MXITST                                               40.41
            CALL MSGERR (1, '[itermx] is replaced; see user manual')      40.41
          ELSE IF (KEYWIS ('NONST')) THEN                                 40.41
            CALL ININTG ('MXITNS' , MXITNS   , 'REQ', 0 )                 40.41
          ENDIF                                                           40.41
          CALL INREAL ('LIMITER', PNUMS(20), 'UNC', 0.)                   40.41
        END IF
!
!       *** numerical scheme in directional space (standard  ***
!       *** option implicit scheme)                          ***
!
        CALL INKEYW ('STA','  ')
        IF (KEYWIS ('DIR')) THEN
          CALL INREAL ('CDD'    , PNUMS(6) , 'UNC', 0.)
          CALL INREAL ('CDLIM'  , PNUMS(17), 'UNC', 0.)                   30.80
          IF (PNUMS(17).LT.0.) IREFR = 1                                  30.80
          IF (PNUMS(17).GT.0.) IREFR = -1                                 40.02
          IF (EQREAL(PNUMS(17),0.)) THEN                                  30.80
            IREFR = 0                                                     30.80
            CALL MSGERR(0, 'Refraction deactivated')                      40.02
          ENDIF
          CALL INKEYW ('STA','    ')                                      41.07
          IF (KEYWIS ('WNUM')) PNUMS(32) = 1.                             41.07
        ENDIF
!
!       limit Ctheta if user want so                                      41.06
!
        CALL INKEYW ('STA','  ')
        IF (KEYWIS ('REFRL')) THEN
          PNUMS(29) = 1.
          CALL INREAL ('FRLIM', PNUMS(26) ,'UNC', 0.)
          CALL INREAL ('POWER', PNUMS(27), 'UNC', 0.)
        ENDIF
!
!       *** numerical scheme in frequency space :                  ***
!       *** 1) Fully implicit scheme in frequency space (iterative ***
!       ***    SIP solver)                                         ***
!       *** 2) Explicit scheme in frequency space:                 ***
!       ***    Energy is removed from the spectrum based on a CFL  ***
!       ***    criterion                                           ***
!       *** 3) Explicit scheme in frequency space:                 ***
!       ***    No CFL limitation near blocking point --> unstable  ***
!       ***    integration. Therefore a filter is applied with a   ***
!       ***    diffusion coeff.                                    ***
!
        CALL INKEYW ('STA','  ')
        IF (KEYWIS('SIGIM') .OR. KEYWIS ('IMP')) THEN                     30.20
!         *** implicit solver ***
!         ***  This is the default option   PNUMS(8) = 1.  ***
          PNUMS(8) = 1.
          CALL INREAL ('CSS'  , PNUMS (7), 'UNC', 0.)
          CALL INREAL ('EPS1' , PNUMS(11), 'UNC', 0.)
          CALL INREAL ('EPS2' , PNUMS(12), 'UNC', 0.)
          CALL INREAL ('OUTP' , PNUMS(13), 'UNC', 0.)
          CALL INREAL ('NITER', PNUMS(14), 'UNC', 0.)
        ELSE IF (KEYWIS('SIGEX') .OR. KEYWIS('EXP')) THEN                 30.20
!
!         *** 2) explicit scheme ***
          PNUMS(8) = 2.
          CALL INREAL ('CSS'  , PNUMS (7), 'UNC', 0.)                     41.07
          CALL INREAL ('CFL'  , PNUMS(19), 'UNC', 0.)
!
        ELSE IF (KEYWIS ('FIL')) THEN
!         *** 3) explicit scheme -> filter the spectrum  ***
          PNUMS(8) = 3.
          CALL INREAL ('DIFFC'  , PNUMS(9), 'UNC', 0.)
!
        END IF
        CALL INKEYW ('STA','  ')                                          30.82
        IF (KEYWIS('SETUP')) THEN                                         30.82
!         *** iterative solver        ***                                 30.82
!         *** Settings for the setup  ***                                 30.82
          CALL INREAL ('EPS2' , PNUMS(23), 'UNC', 0.)                     30.82
          CALL INREAL ('OUTP' , PNUMS(24), 'UNC', 0.)                     30.82
          CALL INREAL ('NITER', PNUMS(25), 'UNC', 0.)                     30.82
        ENDIF                                                             30.82
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!     SETUP     include wave-induced set-up in SWAN calculation
!
! ============================================================
!
!   SETUP [supcor]                                                        30.82
!
! ============================================================
!
      IF (KEYWIS ('SETUP')) THEN                                          32.02
        IF (KSPHER.NE.0) CALL MSGERR (2,
     & 'setup will not be compute correctly with spherical coordinates')  40.13
        IF (PARLL .AND. .NOT.ONED) THEN                                   40.41
           CALL MSGERR(4,'setup is not supported in parallel run')
           RETURN
        END IF
        IF (OPTG.EQ.5) THEN                                               40.80
           CALL MSGERR(2, ' No setup for unstructured grid')              40.80
           LSETUP = 0                                                     40.80
           GOTO 100                                                       40.80
        ENDIF                                                             40.80
        LSETUP = 1                                                        32.02
        CALL INREAL ('SUPCOR', PSETUP(2), 'UNC', 0.)                      30.82
!       *** set pointers for setup and saved depth in array COMPDA ***    32.02
        JSETUP = MCMVAR + 1                                               32.02
        JDPSAV = MCMVAR + 2                                               32.02
        MCMVAR = MCMVAR + 2                                               40.41 32.02
        ALOCMP = .TRUE.                                                   40.97
        GOTO 100                                                          32.02
      ENDIF                                                               32.02
!
!     ------------------------------------------------------------------  40.21
!
!     DIFFRac   include diffraction approximation                         40.21
!                                                                         40.21
! =============================================================           40.21
!                                                                         40.21
!   DIFFRac  [idiffr]  [smpar]  [smnum]  [cgmod]                          40.21
!                                                                         40.21
! =============================================================           40.21
!                                                                         40.21
      IF (KEYWIS ('DIFFR')) THEN                                          40.21
        CALL ININTG ('IDIFFR', IDIFFR, 'STA', 1)                          40.21
        CALL INREAL ('SMPAR', PDIFFR(1), 'STA', 0.0)                      40.21
        CALL INREAL ('SMNUM', PDIFFR(2), 'STA', 0.)                       40.21
        CALL INREAL ('CGMOD', PDIFFR(3), 'STA', 1.)                       40.21
        GOTO 100                                                          40.21
      ENDIF                                                               40.21
!
!     ------------------------------------------------------------------
!
!     SET       setting physical parameters and error counters
!
! =============================================================
!
!   SET  [level]  [nor]  [depmin]  [maxmes]        &
!        [maxerr]  [grav]  [rho] [cdcap] [inrhog]  &
!        [hsrerr]  CARTesian/NAUTical  [pwtail]    &
!        [froudmax]  [printf]  [prtest]
!
! =============================================================
!
      IF (KEYWIS ('SET')) THEN
         CALL INREAL ('LEVEL',  WLEV,   'UNC', 0.)
         CALL INREAL ('NOR',    DNORTH, 'UNC', 0.)
         CALL INREAL ('DEPMIN', DEPMIN, 'UNC', 0.)
         CALL ININTG ('MAXMES', MAXMES, 'UNC', 0)                         30.70
         CALL ININTG ('MAXERR', MAXERR, 'UNC', 0)
         CALL INREAL ('GRAV',   GRAV,   'UNC', 0.)
         CALL INREAL ('RHO',    RHO,    'UNC', 0.)
         CALL INREAL ('CDCAP',  CDCAP,  'UNC', 0.)
         CALL ININTG ('INRHOG', INRHOG, 'UNC', 0)
         IF (INRHOG.EQ.0) THEN                                            30.20
           OVUNIT(7) = 'm2/s'
           OVUNIT(9) = 'm2/s'
           OVUNIT(19) = 'm3/s'
           OVUNIT(21) = 'm2s'
           OVUNIT(22) = 'm2'
           OVUNIT(29) = 'm2'
         ELSE
           OVUNIT(7) = 'W/m2'
           OVUNIT(9) = 'W/m2'
           OVUNIT(19) = 'W/m'
           OVUNIT(21) = 'Js/m2'
           OVUNIT(22) = 'J/m2'
           OVUNIT(29) = 'J/m2'
         ENDIF                                                            10.17
!
         CALL INREAL ('HSRERR', HSRERR, 'UNC', 0.)                        32.01
!
         CALL INKEYW ('STA',' ')                                          32.01
         IF (KEYWIS ('NAUT')) THEN                                        40.02
           BNAUT = .TRUE.                                                 32.01
           OUTPAR(4) = 0.                                                 40.02
         ENDIF                                                            40.02
         IF (KEYWIS ('CART')) BNAUT = .FALSE.                             30.70
         IF ( ITEST .GE. 20 ) THEN                                        32.01
           WRITE (PRTEST,*) ' Set NAUT command: BNAUT=', BNAUT            32.01
         ENDIF                                                            32.01
!
         CALL INREAL ('PWTAIL', PWTAIL(1), 'UNC', 0.)
         IF (CHGVAL) THEN
           IF (PWTAIL(1).LE.1.) CALL MSGERR (3, 'Incorrect PWTAIL')
           PWTAIL(3) = PWTAIL(1) + 1.
         ENDIF
         CALL INREAL ('FROUDMAX', PNUMS(18), 'UNC', 0.)                   30.50
         CALL ININTG ('PRINTF', PRINTF, 'UNC', 0)
         CALL ININTG ('PRTEST', PRTEST, 'UNC', 0)
         GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     QUANTITY  setting parameters for output quantities                  40.03
!
! =============================================================
!
!   QUANTity  <...>  'short'  'long'  [lexp]  [hexp]  [excv]  &           40.03
!
!         [ref]                    {for output quantity TSEC}
!         [power]                  {for output quantity WLEN, PER or RPER}
!         [fswell]                 {for output quantity HSWELL}           40.03
!         [fmin]  [fmax]           {for integral parameters}              40.87
!         PROBLEM/FRAME            {for directions and vectors)           40.03
!
! =============================================================
!
      IF (KEYWIS('QUANT')) THEN                                           40.13
         NVAR    = 0                                                      40.87
         FRSTQ%I = 0                                                      40.87
         NULLIFY(FRSTQ%NEXTI)                                             40.87
         CURRQ => FRSTQ                                                   40.87
 70      CALL SVARTP (IVTYPE)                                             40.87
         IF (IVTYPE.GE.1 .AND. IVTYPE.LE.NMOVAR) THEN                     40.87
            NVAR = NVAR+1                                                 40.87
            ALLOCATE(TMPQ)                                                40.87
            TMPQ%I = IVTYPE                                               40.87
            NULLIFY(TMPQ%NEXTI)                                           40.87
            CURRQ%NEXTI => TMPQ                                           40.87
            CURRQ => TMPQ                                                 40.87
            CALL INCSTR ('SHORT', QOVSNM, 'STA', ' ')                     40.87
            CALL INCSTR ('LONG' , QOVLNM, 'STA', ' ')                     40.87
            CALL INREAL ('LEXP' , QR(1) , 'STA',  999. )                  40.87
            CALL INREAL ('HEXP' , QR(2) , 'STA', -999. )                  40.87
            CALL INREAL ('EXCV' , QR(3) , 'STA',  999. )                  40.87
            CALL INCTIM (ITMOPT, 'REF', QR(4), 'STA', -999.)              40.87
            CALL INREAL ('POWER', QR(5) , 'STA', -999.)                   40.87
            CALL INREAL ('FSWELL', QR(6), 'STA', -999.)                   40.87
            CALL INREAL ('FMIN' , QR(8) , 'STA',  999. )                  40.87
            CALL INREAL ('FMAX' , QR(9) , 'STA', -999. )                  40.87
            CALL INKEYW ('STA', ' ')                                      40.87
            IF (KEYWIS('PROBLEM') .OR. KEYWIS('USER')) THEN               40.87
               QR(7) = 0.                                                 40.87
            ELSE IF (KEYWIS('FRAME')) THEN                                40.87
               QR(7) = 1.                                                 40.87
            ELSE                                                          40.87
               QR(7) = -999.                                              40.87
            ENDIF                                                         40.87
            GOTO 70                                                       40.87
         ELSEIF (IVTYPE.NE.99) THEN                                       40.87
           CALL MSGERR (2, 'unknown quantity ')                           40.03
         ENDIF                                                            40.87
!
         IF (NVAR.GT.0) THEN                                              40.87
            CURRQ => FRSTQ%NEXTI                                          40.87
            DO JJ = 1, NVAR                                               40.87
               IVTYPE = CURRQ%I                                           40.87
               IF (QOVSNM.NE.' ') OVSNAM(IVTYPE) = QOVSNM                 40.87
               IF (QOVLNM.NE.' ') OVLNAM(IVTYPE) = QOVLNM                 40.87
               IF (QR(1).NE. 999.) OVLEXP(IVTYPE) = QR(1)                 40.87
               IF (QR(2).NE.-999.) OVHEXP(IVTYPE) = QR(2)                 40.87
               IF (QR(3).NE. 999.) OVEXCV(IVTYPE) = QR(3)                 40.87
               IF (IVTYPE.EQ.40 .OR. IVTYPE.EQ.41) THEN
                  IF (NSTATM.EQ.0) CALL MSGERR (2,
     &                           'Time output asked in stationary mode')
               ENDIF
               IF (IVTYPE.EQ.41) THEN
                  IF (QR(4).NE.-999.) OUTPAR(1) = QR(4)                   40.87
               ELSE IF (IVTYPE.EQ.42 .OR. IVTYPE.EQ.43) THEN
                  IF (QR(5).NE.-999.) OUTPAR(2) = QR(5)                   40.87
               ELSE IF (IVTYPE.EQ.17) THEN
                  IF (QR(5).NE.-999.) OUTPAR(3) = QR(5)                   40.87
               ELSE IF (IVTYPE.EQ.44) THEN
                  IF (QR(6).NE.-999.) OUTPAR(5) = QR(6)                   40.87
               ENDIF
               IF ( ANY( IVTYPE == IVOTP ) ) THEN                         40.87
                  INDX = MINLOC(IVOTP, IVOTP==IVTYPE)                     40.87
                  IF (QR(8).NE.999.) THEN                                 40.87
                     OUTPAR(INDX(1)+ 5) = 1.                              40.87
                     OUTPAR(INDX(1)+20) = QR(8)                           40.87
                  ENDIF                                                   40.87
                  IF (QR(9).NE.-999.) THEN                                40.87
                     OUTPAR(INDX(1)+ 5) = 1.                              40.87
                     OUTPAR(INDX(1)+35) = QR(9)                           40.87
                  ENDIF                                                   40.87
                  IF ( OUTPAR(INDX(1)+20).GT.OUTPAR(INDX(1)+35) )
     &               CALL MSGERR (2,'fmin cannot be larger than fmax')    40.87
               ELSE                                                       40.87
                  IF (QR(8).NE.999. .OR. QR(9).NE.-999.) CALL MSGERR (1,  40.87
     &                'Integration range not allowed for this quantity')  40.87
               ENDIF                                                      40.87
               IF (OVSVTY(IVTYPE).EQ.2 .OR. OVSVTY(IVTYPE).EQ.3) THEN
!              direction or vector
                  IF (QR(7).NE.-999.) THEN
                     OUTPAR(4) = QR(7)                                    40.87
                     IF (BNAUT.AND.OVSVTY(IVTYPE).EQ.2) CALL MSGERR (1,
     &                    'option not allowed with Nautical convention')
                  ENDIF
               ENDIF
               CURRQ => CURRQ%NEXTI                                       40.87
            END DO
            DEALLOCATE(TMPQ)                                              40.87
         ENDIF
         GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     BREAK     parameters surf breaking
!
! ============================================================
!
!           | -> CONstant [alpha] [gamma]                                      |
!           |                                                                  |
!           |    VARiable [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
! BREaking <                                                                    > &
!           |    RUEssink [alpha] [a] [b]                                      |
!           |                                                                  |
!           |    TG       [alpha] [gamma] [pown]                               |
!
!
!       ( FREQDep [power] [fmin] [fmax] ) (fmin and fmax not documented)
!
! ============================================================
!
      IF (KEYWIS ('BRE')) THEN
        CALL INKEYW ('STA', 'CON')
        IF (KEYWIS('CON')) THEN
          ISURF = 1
          CALL INREAL ('ALPHA', PSURF(1), 'STA', 1.0)
          CALL INREAL ('GAMMA', PSURF(2), 'STA', 0.73)
        ELSE IF (KEYWIS('VAR') .OR. KEYWIS('NEL')) THEN                   41.03
          ISURF = 2
          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.5)
          CALL INREAL ('GAMMIN', PSURF(4), 'STA', 0.55)
          CALL INREAL ('GAMMAX', PSURF(5), 'STA', 0.81)
          CALL INREAL ('GAMNEG', PSURF(6), 'STA', 0.73)
          CALL INREAL ('COEFF1', PSURF(7), 'STA', 0.88)
          CALL INREAL ('COEFF2', PSURF(8), 'STA', 0.012)
        ELSE IF (KEYWIS('RUE')) THEN                                      41.03
          ISURF = 3
          CALL INREAL ('ALPHA', PSURF(1), 'STA', 1.0)
          CALL INREAL ('A'    , PSURF(4), 'STA', 0.76)
          CALL INREAL ('B'    , PSURF(5), 'STA', 0.29)
        ELSE IF (KEYWIS('TG')) THEN                                       41.03
          ISURF = 4
          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.0)
          CALL INREAL ('GAMMA',  PSURF(2), 'STA', 0.42)
          CALL INREAL ('POWN' ,  PSURF(4), 'STA', 4.0)
        ENDIF
!
        CALL INKEYW ('STA', '  ')                                         41.06
        IF (KEYWIS ('FREQD')) THEN                                        41.06
           IFRSRF = 1                                                     41.06
           CALL INREAL ('POWER', PSURF(11), 'STA', 2.0)                   41.06
           CALL INREAL ('FMIN' , PSURF(12), 'STA', 0.0)                   41.06
           CALL INREAL ('FMAX' , PSURF(13), 'STA', 1000.)                 41.06
        END IF                                                            41.06
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     WCAP        parameters whitecapping (NOT documented)
!
! =============================================================
!
!        | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]                34.00
!        |
!        |    JANSsen [cds1]  [delta] [pwtail]
!        |
!        |    LHIG    [cflhig]
!        |
! WCAP  <     BJ      [bjstp] [bjalf]
!        |
!        |    KBJ     [bjstp] [bjalf] [kconv]
!        |                                                                40.30
!        |    AB      [cds2] [br] [p0] [powst] [powk]                     40.53
!
! =============================================================
!
      IF (KEYWIS ('WCAP')) THEN
        CALL INKEYW ('STA','KOM')                                         970220
        IF (KEYWIS ('KOM')) THEN
!         *** whitecapping according to Komen et al. (1984) ***
          IWCAP = 1
          CALL INREAL ('CDS2',  PWCAP(1),  'UNC', 0.)                     20.73
          CALL INREAL ('STPM',  PWCAP(2),  'UNC', 0.)                     20.73
          CALL INREAL ('POWST', PWCAP(9),  'UNC', 0.)                     34.00
          CALL INREAL ('DELTA', PWCAP(10), 'UNC', 0.)                     34.00
          CALL INREAL ('POWK',  PWCAP(11), 'UNC', 0.)                     34.00
        ELSE IF ( KEYWIS ('JANS')) THEN
!         *** whitecapping according to Janssen (1989, 1991) ***
          IWCAP = 2
          CALL INREAL ('CDS1',  PWCAP(3), 'UNC', 0.)                      20.73
          CALL INREAL ('DELTA', PWCAP(4), 'UNC', 0.)
!
!         Recalculate coefficients that are actually used in the          40.02
!         whitecapping routines                                           40.02
          PWCAP(1)  = PWCAP(3) * (PWCAP(2) ** PWCAP(9))                   40.02
          PWCAP(10) = PWCAP(4)                                            40.02
          JUSTAR = MCMVAR+1                                               30.22
          JZEL   = MCMVAR+2                                               30.22
          JCDRAG = MCMVAR+3                                               30.22
          JTAUW  = MCMVAR+4                                               30.22
          MCMVAR = MCMVAR+4                                               30.22
          ALOCMP = .TRUE.                                                 40.97
          CALL INREAL ('PWTAIL', PWTAIL(1), 'STA', 5.)                    30.70
          IF (PWTAIL(1).LE.1.) CALL MSGERR (3, 'Incorrect PWTAIL')        30.70
          PWTAIL(3) = PWTAIL(1) + 1.                                      30.70
        ELSE IF ( KEYWIS ('LHIG')) THEN
!         *** whitecapping according to Longuett Higgins ***
          IWCAP = 3
          CALL INREAL ('CFLHIG', PWCAP(5), 'UNC', 0.)
        ELSE IF ( KEYWIS ('BJ')) THEN
!         *** whitecapping according to Battjes/Janssen formulation ***
          IWCAP = 4
          CALL INREAL ('BJSTP' , PWCAP(6), 'UNC', 0.)
          CALL INREAL ('BJALF' , PWCAP(7), 'UNC', 0.)
        ELSE IF ( KEYWIS ('KBJ')) THEN
!         *** whitecapping according to a combination of Komen et al ***
!         *** and Battjes/Janssen                                    ***
          IWCAP = 5
          CALL INREAL ('BJSTP' , PWCAP(6), 'UNC', 0.)
          CALL INREAL ('BJALF' , PWCAP(7), 'UNC', 0.)
          CALL INREAL ('KCONV' , PWCAP(8), 'UNC', 0.)
        ELSE IF ( KEYWIS ('AB')) THEN                                     40.53
!         *** whitecapping according to Alves and Banner (2003) ***       40.53
          IWCAP = 7                                                       40.53
          CALL INREAL ('CDS2',  PWCAP(1),  'STA', 5.0E-5)                 40.53
          CALL INREAL ('BR',    PWCAP(12), 'STA', 1.75E-3)                40.53
          CALL INREAL ('P0',    PWCAP(10), 'STA', 4.)                     40.53
          CALL INREAL ('POWST', PWCAP(9),  'STA', 0.)                     40.53
          CALL INREAL ('POWK',  PWCAP(11), 'STA', 0.)                     40.53
        ELSE
          CALL WRNKEY
        ENDIF
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     FRIC  setting parameters for bottom friction
!
! ===============================================
!
!                   |              | -> CONstant [cfjon]
!                   | -> JONswap  <
!                   |              |    VARiable [cfj1] [cfj2] [dsp1] [dsp2]   (NOT documented)
!                   |
!                   |    COLLins   [cfw]    &
!   FRICtion       <
!                   |              [cfc]      (NOT documented)
!                   |
!                   |    MADsen    [kn]
!
! ===============================================
!
      IF (KEYWIS ('FRIC')) THEN
        IBOT = 1
        CALL INKEYW ('STA','JON')
        IF (KEYWIS('JON')) THEN
          CALL INKEYW ('STA', 'CON')
          IF (KEYWIS('VAR')) THEN                                         41.04
            IBOT = 4
            CALL INREAL ('CFJ1', PBOT(6), 'STA', 0.038)
            CALL INREAL ('CFJ2', PBOT(7), 'STA', 0.067)
            CALL INREAL ('DSP1', PBOT(8), 'STA', 10.)
            CALL INREAL ('DSP2', PBOT(9), 'STA', 30.)
          ELSE
            CALL IGNORE ('CON')
            IBOT = 1
            CALL INREAL('CFJON',PBOT(3),'UNC',0.)
          ENDIF
        ELSE IF (KEYWIS('COLL')) THEN                                     20.68
          IBOT = 2
          CALL INREAL('CFW',PBOT(2),'UNC',0.)
          CALL INREAL('CFC',PBOT(1),'UNC',0.)
        ELSE IF (KEYWIS('MAD')) THEN
          IBOT = 3
          CALL INREAL('KN',PBOT(5),'UNC',0.)
        ELSE
          CALL WRNKEY
        ENDIF
        GOTO 100
      ENDIF
!
!   ------------------------------------------------------------------
!   VEGEtation   < [height]  [diamtr]  [nstems]  [drag] >
!   ------------------------------------------------------------------
      IF (KEYWIS ('VEGE')) THEN
        IVEG  = 1
        ILMAX = 0
        FRSTV%H = 0.
        FRSTV%D = 0.
        FRSTV%N = 0
        FRSTV%C = 0.
        NULLIFY(FRSTV%NEXTV)
        CURRV => FRSTV
 303    CALL INREAL ('HEIGHT',HH,'REP',-1.)
        IF ( HH.NE.-1. ) THEN
           IF ( HH.LT.0. ) THEN
              CALL MSGERR (2,'height is negative')
              HH = 0.
           END IF
           CALL INREAL ('DIAMTR',DD,'REQ', 0.)
           IF ( DD.LT.0. ) THEN
              CALL MSGERR (2,'stem diameter is negative')
              DD = 0.
           END IF
           CALL ININTG ('NSTEMS',NN,'REQ', 1 )
           IF ( NN.LE.0 ) THEN
              CALL MSGERR (2,'number of stems is negative or zero')
              NN = 1
           END IF
           CALL INREAL ('DRAG'  ,CC,'REQ', 0.)
           IF ( CC.LT.0. ) THEN
              CALL MSGERR (2,'drag coefficient is negative')
              CC = 0.
           END IF
           ILMAX = ILMAX+1
           ALLOCATE(TMPV)
           TMPV%H = HH
           TMPV%D = DD
           TMPV%N = NN
           TMPV%C = CC
           NULLIFY(TMPV%NEXTV)
           CURRV%NEXTV => TMPV
           CURRV => TMPV
           GO TO 303
        END IF
        IF (ILMAX.EQ.0) CALL MSGERR(2,'No vegetation parameters found')
        IF (.NOT.ALLOCATED(LAYH  )) ALLOCATE(LAYH  (ILMAX))
        IF (.NOT.ALLOCATED(VEGDIL)) ALLOCATE(VEGDIL(ILMAX))
        IF (.NOT.ALLOCATED(VEGNSL)) ALLOCATE(VEGNSL(ILMAX))
        IF (.NOT.ALLOCATED(VEGDRL)) ALLOCATE(VEGDRL(ILMAX))
        CURRV => FRSTV%NEXTV
        DO JJ = 1, ILMAX
           LAYH  (JJ) = CURRV%H
           VEGDIL(JJ) = CURRV%D
           VEGNSL(JJ) = REAL(CURRV%N)
           VEGDRL(JJ) = CURRV%C
           CURRV => CURRV%NEXTV
        END DO
        DEALLOCATE(TMPV)
        GOTO 100
      END IF
!
!     ------------------------------------------------------------------
!
!     WIND      parameters uniform wind field
!
! ==========================================
!
!   WIND  [vel]  [dir]  [astd]                                            40.03
!
! ==========================================
!
      IF (KEYWIS ('WIND')) THEN
        LWINDR = 1                                                        30.10
        IWIND  = LWINDM                                                   30.10
        VARWI = .FALSE.
        CALL INREAL('VEL',U10,'REQ',0.)
        CALL INREAL('DIR',WDIP,'REQ',0.)
        CALL INREAL('ASTD',CASTD,'STA',0.)                                40.03
!
!       *** Convert (if necessary) WDIP from nautical degrees ***         32.01
!       *** to cartesian degrees                              ***         32.01
!
        WDIP = DEGCNV (WDIP)                                              32.01
!
        IF (IWIND.EQ.0) IWIND = 4                                         20.74
        ALTMP = WDIP / 360.
        WDIP = PI2 * (ALTMP - NINT(ALTMP))
!
        GOTO 100
      ENDIF
!     ------------------------------------------------------------------
!
!     GEN1        deep water with 1st generation                          20.86
!
! ==========================================
!
! GEN1     [cf10]  [cf20]  [cf30]  [cf40]  [edmlpm]  [cdrag]  [umin]  [cfpm]
!
! ==========================================
!
      IF (KEYWIS('GEN1')) THEN
!       *** initialize first generation model ***
        IGEN = 1                                                          32.06
        LWINDM = 1                                                        30.10
        IF (LWINDR .GT. 0) IWIND = LWINDM                                 30.10
        IQUAD = 0
        IWCAP = 0
!       *** set value for the windparameters ***
!       *** first generation wind model ***
        CALL INREAL ('CF10', PWIND(1), 'UNC', 0.)
        CALL INREAL ('CF20', PWIND(2), 'UNC', 0.)
        CALL INREAL ('CF30', PWIND(3),'UNC',0.)
        CALL INREAL ('CF40', PWIND(4),'UNC',0.)
        CALL INREAL ('EDMLPM', PWIND(10), 'UNC', 0.)
        CALL INREAL ('CDRAG',  PWIND(11), 'UNC', 0.)
        CALL INREAL ('UMIN', PWIND(12), 'UNC', 0.)
        CALL INREAL ('CFPM',  PWIND(13), 'UNC', 0.)
!       if [pwtail] is not changed, make it 5
        IF (EQREAL(PWTAIL(1),4.)) THEN                                    40.00
          PWTAIL(1) = 5.                                                  40.00
          PWTAIL(3) = PWTAIL(1) + 1.                                      30.70
        ENDIF
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     GEN2        deep water with 2nd generation                          20.86
!
! ==========================================
!
! GEN2 [cf10] [cf20] [cf30] [cf40] [cf50] [cf60] [edmlpm] [cdrag] [umin] [cfpm]
!
! ==========================================
!
      IF (KEYWIS('GEN2')) THEN
!       *** initialize second generation model ***
        IGEN = 2                                                          32.06
        LWINDM = 2                                                        30.10
        IF (LWINDR .GT. 0) IWIND = LWINDM                                 30.10
        IQUAD = 0
        IWCAP = 0
!       *** set value for the windparameters ***
        CALL INREAL ('CF10', PWIND(1), 'UNC', 0.)
        CALL INREAL ('CF20', PWIND(2), 'UNC', 0.)
        CALL INREAL ('CF30', PWIND(3),'UNC',0.)
        CALL INREAL ('CF40', PWIND(4),'UNC',0.)
        CALL INREAL ('CF50', PWIND(5),'UNC',0.)
        CALL INREAL ('CF60', PWIND(6), 'UNC', 0.)
        CALL INREAL ('EDMLPM', PWIND(10), 'UNC', 0.)
        CALL INREAL ('CDRAG',  PWIND(11), 'UNC', 0.)
        CALL INREAL ('UMIN', PWIND(12), 'UNC', 0.)
        CALL INREAL ('CFPM',  PWIND(13), 'UNC', 0.)
!       if [pwtail] is not changed, make it 5
        IF (EQREAL(PWTAIL(1),4.)) THEN                                    40.00
          PWTAIL(1) = 5.                                                  40.00
          PWTAIL(3) = PWTAIL(1) + 1.                                      30.70
        ENDIF
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     GEN3        deep water with 3d generation                           20.86
!
! ======================================================================================
!
!       |    JANSsen [cds1] [delta]                        |              40.02
!       |                                                  |
!       | -> KOMen   [cds2] [stpm]                         |              34.00
! GEN3 <                                                    > (AGROW [a])
!       |    YAN     (NOT documented)                      |
!       |                                                  |
!       |    WESTHuysen [cds2] [br] [p0] [powst] [powk]    |              40.53
!
! ======================================================================================
!
      IF (KEYWIS('GEN3')) THEN
!       *** initialize third generation model ***
        IGEN = 3                                                          32.06
        CALL INKEYW ('STA', 'KOM')                                        970220
        IF (KEYWIS('JANS')) THEN                                          970220
          LWINDM = 4                                                      30.1x
          IF (LWINDR .GT. 0) IWIND = LWINDM                               30.1x
!         *** whitecapping according to Janssen (1989, 1991) ***
          IWCAP = 2
          CALL INREAL ('CDS1',  PWCAP(3), 'UNC', 0.)                      20.73
          CALL INREAL ('DELTA', PWCAP(4), 'UNC', 0.)
!
!         Recalculate coefficientes that are actually used in the         40.02
!         whitecapping routines                                           40.02
          PWCAP(1)  = PWCAP(3) * (PWCAP(2) ** PWCAP(9))                   40.02
          PWCAP(10) = PWCAP(4)                                            40.02
          JUSTAR = MCMVAR+1                                               30.72
          JZEL   = MCMVAR+2                                               30.72
          JCDRAG = MCMVAR+3                                               30.72
          JTAUW  = MCMVAR+4                                               30.72
          MCMVAR = MCMVAR+4                                               30.72
          ALOCMP = .TRUE.                                                 40.97
!         if [pwtail] is not changed, make it 5
          IF (EQREAL(PWTAIL(1),4.)) THEN                                  40.00
            PWTAIL(1) = 5.                                                40.00
            PWTAIL(3) = PWTAIL(1) + 1.                                    30.70
          ENDIF
        ELSE IF (KEYWIS('YAN')) THEN
!         option not documented in user manual
          LWINDM = 5                                                      30.1x
          IF (LWINDR .GT. 0) IWIND = LWINDM                               30.1x
        ELSE IF (KEYWIS('WESTH')) THEN                                    40.53
!         *** wind according to (adapted) Yan (1987) ***                  40.53
          LWINDM = 5                                                      40.53
          IF (LWINDR .GT. 0) IWIND = LWINDM                               40.53
!         *** whitecapping according to Alves and Banner (2003) ***       40.53
          IWCAP = 7                                                       40.53
          CALL INREAL ('CDS2',  PWCAP(1),  'STA', 5.0E-5)                 40.53
          CALL INREAL ('BR',    PWCAP(12), 'STA', 1.75E-3)                40.53
          CALL INREAL ('P0',    PWCAP(10), 'STA', 4.)                     40.53
          CALL INREAL ('POWST', PWCAP(9),  'STA', 0.)                     40.53
          CALL INREAL ('POWK',  PWCAP(11), 'STA', 0.)                     40.53
        ELSE
          CALL IGNORE ('KOM')                                             970220
          LWINDM = 3                                                      30.1x
          IF (LWINDR .GT. 0) IWIND = LWINDM                               30.1x
!         *** whitecapping according to Komen et al. (1984) ***
          IWCAP = 1
          CALL INREAL ('CDS2', PWCAP(1), 'UNC', 0.)                       20.73
          CALL INREAL ('STPM', PWCAP(2), 'UNC', 0.)                       20.73
        ENDIF
        CALL INKEYW ('STA', ' ')
        IF (KEYWIS('AGROW')) THEN                                         7/MAR
          CALL INREAL ('A',PWIND(31),'STA',0.0015)
        ELSE                                                              30.60
          IF (NSTATM.EQ.1 .AND. ICOND.EQ.0) ICOND = 1                     30.70
        ENDIF
        GOTO 100
      ENDIF
!     ----------------------------------------------------------------    40.17

!                               | CNL4 < [cnl4] >               |         40.17
!     MDIA LAMbda < [lambda] > <                                 >        40.17
!                               | CNL4_12 < [cnl4_1] [cnl4_2] > |         40.17

!     ----------------------------------------------------------------    40.17

      IF (KEYWIS('MDIA')) THEN                                            40.17
        IF (ALLOCATED(LAMBDA)) THEN                                       40.17
          DEALLOCATE (LAMBDA,CNL4_1,CNL4_2)                               40.17
        ENDIF                                                             40.17
        ALLOCATE (RLAMBDA(1000))                                          40.17
        CALL INKEYW ('STA', '   ')                                        40.17
        IF (KEYWIS('LAM')) THEN                                           40.17
          MORE = .TRUE.                                                   40.17
          ILAMBDA = 0                                                     40.17
          DO WHILE (MORE)                                                 40.17
            CALL INREAL('LAMBDA',RLAMBDA(ILAMBDA+1),'STA',-1.)            40.17
            IF (RLAMBDA(ILAMBDA+1).LT.0.) THEN                            40.17
              MORE = .FALSE.                                              40.17
            ELSE                                                          40.17
              ILAMBDA = ILAMBDA + 1                                       40.17
            ENDIF                                                         40.17
          ENDDO                                                           40.17
          MDIA = ILAMBDA                                                  40.17
          ALLOCATE (LAMBDA(MDIA))                                         40.17
          LAMBDA(1:MDIA) = RLAMBDA(1:MDIA)                                40.17
          DEALLOCATE (RLAMBDA)                                            40.17
        ELSE                                                              40.17
          CALL WRNKEY                                                     40.17
        END IF                                                            40.17
        CALL INKEYW ('STA', '   ')                                        40.17
        ALLOCATE (CNL4_1(MDIA),CNL4_2(MDIA))                              40.17
        IF (KEYWIS('CNL4_12')) THEN                                       40.17
          DO ICNL4=1,MDIA                                                 40.17
            CALL INREAL('CNL4_1',CNL4_1(ICNL4),'REQ',0.)                  40.17
            CALL INREAL('CNL4_2',CNL4_2(ICNL4),'REQ',0.)                  40.17
          END DO                                                          40.17
        ELSEIF (KEYWIS('CNL4')) THEN                                      40.17
          DO ICNL4=1,MDIA                                                 40.17
            CALL INREAL('CNL4',CNL4_1(ICNL4),'REQ',0.)                    40.17
          END DO                                                          40.17
          CNL4_2 = CNL4_1                                                 40.17
        ELSE                                                              40.17
          CALL WRNKEY                                                     40.17
        END IF                                                            40.17
        CNL4_1 = CNL4_1 * ((2.*PI)**9)                                    40.17
        CNL4_2 = CNL4_2 * ((2.*PI)**9)                                    40.17
        GOTO 100                                                          40.17
      ENDIF                                                               40.17

!     ------------------------------------------------------------------
!
!     GROWTH       parameters wind input source term                      20.67
!     Note: command does not appear in user manual for SWAN               20.86
!
! ==========================================
!
!         |   G1 [cf10] [cf20] [cf30] [cf40]               & |
!         |                                                  |
!         |      [edmlpm] [cdrag] [umin] [cfpm]              |
!         |                                                  |
!         |   G2 [cf10] [cf20] [cf30] [cf40] [cf50] [cf60] & |
!         |                                                  |
! GROWTH <       [edmlpm] [cdrag] [umin] [cfpm]              | (NOT documented)
!         |                                                  |
!         |       |   JANSsen [pwtail] |                     |
!         | ->G3 <                      > (AGROW [a])        |
!         |       | ->KOMen            |                     |
!         |       |                    |                     |
!         |       |   YAN              |                     |
!
!
      IF (KEYWIS('GROWTH') .OR. KEYWIS ('PARW')) THEN                     20.67
!         *** initialize first generation model ***
          CALL INKEYW ('STA', 'G3')                                       20.6x
          IF (KEYWIS('G1') .OR. KEYWIS ('SNY1')) THEN                     20.6x
            IGEN = 1                                                      32.06
            LWINDM = 1                                                    40.13
            IQUAD = 0                                                     20.RR
!           *** set value for the windparameters ***
!           *** first generation wind model ***
            CALL INREAL ('CF10', PWIND(1), 'UNC', 0.)
            CALL INREAL ('CF20', PWIND(2), 'UNC', 0.)
            CALL INREAL ('CF30', PWIND(3),'UNC',0.)
            CALL INREAL ('CF40', PWIND(4),'UNC',0.)
            CALL INREAL ('CF50', PWIND(5),'UNC',0.)
            CALL INREAL ('CF60', PWIND(6), 'UNC', 0.)
            CALL INREAL ('CF70', PWIND(7), 'UNC', 0.)
            CALL INREAL ('CF80', PWIND(8), 'UNC', 0.)
            CALL INREAL ('RHOAW', PWIND(9), 'UNC', 0.)
            CALL INREAL ('EDMLPM', PWIND(10), 'UNC', 0.)
            CALL INREAL ('CDRAG',  PWIND(11), 'UNC', 0.)
            CALL INREAL ('UMIN', PWIND(12), 'UNC', 0.)
            CALL INREAL ('CFPM',  PWIND(13), 'UNC', 0.)
          ELSE IF (KEYWIS('G2') .OR. KEYWIS('SNY2')) THEN                 20.6x
            IGEN = 2                                                      32.06
            LWINDM = 2                                                    40.13
            IQUAD = 0                                                     20.RR
!           *** second generation wind model ***
            CALL INREAL ('CF10', PWIND(1), 'UNC', 0.)
            CALL INREAL ('CF20', PWIND(2), 'UNC', 0.)
            CALL INREAL ('CF30', PWIND(3),'UNC',0.)
            CALL INREAL ('CF40', PWIND(4),'UNC',0.)
            CALL INREAL ('CF50', PWIND(5),'UNC',0.)
            CALL INREAL ('CF60', PWIND(6), 'UNC', 0.)
            CALL INREAL ('CF70', PWIND(7), 'UNC', 0.)
            CALL INREAL ('CF80', PWIND(8), 'UNC', 0.)
            CALL INREAL ('RHOAW', PWIND(9), 'UNC', 0.)
            CALL INREAL ('EDMLPM', PWIND(10), 'UNC', 0.)
            CALL INREAL ('CDRAG',  PWIND(11), 'UNC', 0.)
            CALL INREAL ('UMIN', PWIND(12), 'UNC', 0.)
            CALL INREAL ('CFPM',  PWIND(13), 'UNC', 0.)
          ELSE IF (KEYWIS('G3')) THEN
            IGEN = 3                                                      32.06
            CALL INKEYW ('STA', 'KOM')
            IF (KEYWIS('JANS')) THEN
              LWINDM = 4                                                  40.13
              CALL INREAL ('PWTAIL', PWTAIL(1), 'STA', 5.)                30.70
              IF (PWTAIL(1).LE.1.) CALL MSGERR (3, 'Incorrect PWTAIL')    30.70
              PWTAIL(3) = PWTAIL(1) + 1.                                  30.70
            ELSE IF (KEYWIS('YAN')) THEN
              LWINDM = 5                                                  40.13
            ELSE
              CALL IGNORE ('KOM')
              LWINDM = 3                                                  40.13
            ENDIF
            CALL INKEYW ('STA', ' ')
            IF (KEYWIS('AGROW')) THEN                                     7/MAR
              CALL INREAL ('A',PWIND(31),'STA',0.0015)
            ELSE                                                          30.60
              IF (NSTATM.EQ.1 .AND. ICOND.EQ.0) ICOND = 1                      30.70
            ENDIF
          ELSE
            CALL WRNKEY
          END IF
          IF (LWINDR.GT.0) IWIND = LWINDM                                 40.13
          GOTO 100
      ENDIF
!     ------------------------------------------------------------------
!     ***  parameters nonlinear 4 wave interactions ***
      IF (KEYWIS ('QUAD')) THEN
!
! ===================================================================
!
!  QUADrupl [iquad] [lambda] [cnl4] [csh1] [csh2] [csh3]                  34.00
!
! ===================================================================
!
        CALL ININTG ('IQUAD', IQUAD, 'UNC',  0)
        CALL INREAL ('LAMBDA', PQUAD(1), 'UNC', 0.0)                      34.00
        CALL INREAL ('CNL4', PQUAD(2), 'UNC', 0.0)                        34.00
        CALL INREAL ('CSH1', PQUAD(3), 'UNC', 0.0)                        34.00
        CALL INREAL ('CSH2', PQUAD(4), 'UNC', 0.0)                        34.00
        CALL INREAL ('CSH3', PQUAD(5), 'UNC', 0.0)                        34.00
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     ***  parameters nonlinear 3 wave interactions
!
! ============================================================
!
!  TRIad [trfac] [cutfr] [urcrit] [urslim]                                40.56
!
! ============================================================
!
      IF (KEYWIS ('TRI')) THEN
        ITRIAD = 1                                                        20.81
        CALL INREAL ('TRFAC', PTRIAD(1), 'UNC', 0.0)
        CALL INREAL ('CUTFR', PTRIAD(2), 'STA', 2.5)                      40.61 40.56 30.82
        CALL INREAL ('URCRIT', PTRIAD(4), 'UNC', 0.0)                     40.13
        CALL INREAL ('URSLIM', PTRIAD(5), 'STA', 0.01)                    40.41 40.23 40.13
        GOTO 100
      ENDIF
!
!     ------------------------------------------------------------------
!
!     OFF      switching standard options off
!
! ============================================================
!
!        |  BREaking
!        |
!        |  WCAPping
!        |
!        |  REFrac
! OFF   <
!        |  FSHift
!        |
!        |  QUADrupl
!        |
!        |  WINDGrowth                                                    30.70
!        |
!        |  BNDCHK                                                        40.00
!        |
!        |  RESCALE        (NOT documented)                               40.80
!        |
!        |  SOURCES        (NOT documented)                               40.80
!
! ============================================================
!
      IF (KEYWIS ('OFF')) THEN
        CALL INKEYW ('REQ', ' ')
        IF (KEYWIS ('REF')) THEN
          IREFR = 0
        ELSE IF (KEYWIS ('FSH')) THEN                                     30.20
          ITFRE = 0
        ELSE IF (KEYWIS ('BRE')) THEN
          ISURF = 0
        ELSE IF (KEYWIS ('WCAP')) THEN
          IWCAP = 0
        ELSE IF (KEYWIS ('QUAD')) THEN                                    20.86
          IQUAD = 0                                                       20.86
        ELSE IF (KEYWIS ('WINDG')) THEN                                   30.70
          IWIND = 0                                                       30.70
        ELSE IF (KEYWIS ('BNDCHK')) THEN
!         switch off checking of Hs on boundary                           40.00
          BNDCHK = .FALSE.                                                40.00
        ELSE IF (KEYWIS ('RESCALE')) THEN                                 40.00
!         switch off rescaling solution                                   40.00
          BRESCL = .FALSE.                                                40.00
        ELSE IF (KEYWIS ('SOURCES')) THEN                                 40.80
          OFFSRC = .TRUE.                                                 40.80
          IWIND  = 0                                                      40.80
          IQUAD  = 0                                                      40.80
          IWCAP  = 0                                                      40.80
          ISURF  = 0                                                      40.80
          ITRIAD = 0                                                      40.80
          IBOT   = 0                                                      40.80
        ELSE
          CALL WRNKEY                                                     30.60
        ENDIF
        GOTO 100
      ENDIF
!     ------------------------------------------------------------------
!
!     In case of an empty line in the command file proceed to next line
!
      IF (KEYWRD .EQ. '    ') GOTO 100                                    40.00
!
!     process output requests                                             40.00
!
      CALL SPROUT(FOUND, DEPTH, WLEVL)                                    40.31
      IF (STPNOW()) RETURN                                                34.01
      IF (.NOT.FOUND) THEN
        LEVERR = MAX(LEVERR,3)
        CALL WRNKEY
      ENDIF
      GOTO 100
!
!   * end of subroutine SWREAD *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SINPGR (IGRID1, IGRID2, SNAMEG)                          40.31
!                                                                      *
!***********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.31
      USE SwanGriddata                                                    40.80
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
!     30.60, 30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     32.02: Roeland Ris & Cor van der Schelde
!     30.75, 31.04, 40.13: Nico Booij
!     34.01: Jeroen Adema
!     40.02: IJsbrand Haagsma
!     40.12: IJsbrand Haagsma
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     00.00, Jan. 92
!     30.60, July 97: exception value introduced
!     30.60, July 97: default values for STAGRX, STAGRY, MXINP, MYINP
!     30.60, Aug. 97: format 6021 changed
!     30.60, Aug. 97: keyword EXC not required
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.72, Nov. 97: Header renewed, Common Blocks updated
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     32.02, Feb. 98: Introduced 1D-version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: warning concerning DYINP removed; default MYINP changed
!     30.75, Mar. 98: correction nonstationary wind and current field
!     30.80, Apr. 98: default value for MYINP (curvil. grid) was lost
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Oct. 00: Avoided real/int conflict by introducing replacing
!                     RPOOL for POOL in DPPUTR
!     40.12, Feb. 01: Avoided type conflict for OUTPS
!     40.13, Aug. 01: [xpinp] and [ypinp] are now required in case of
!                     spherical coordinates
!                     Arguments XCGRID and YCGRID removed (not used)
!     40.31, Dec. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Dec. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Read parameters of an input grid
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IGRID1:   input  int grid number for which parameters are read
!     IGRID2:   input  int input grid number for which parameters are read
!                          only relevant if >0                            10.26
!     SNAMEG:   input char name of output frame corresponding to input grid
!
!  9. Subroutines calling
!
!     SWREAD:  Reading and processing of the user commands describing the model
!
!  8. Subroutines used
!
      LOGICAL STPNOW                                                      34.01
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!       subroutine computes the following common data:
!       XCGMIN, XCGMAX, YCGMIN, YCGMAX                                    40.00
!
! 12. Structure
!
!     ---
!
! 13. SOURCE TEXT
!
      INTEGER   IGRID1
      CHARACTER SNAMEG *8
      LOGICAL   KEYWIS                                                    30.00
      TYPE(OPSDAT), POINTER :: OPSTMP                                     40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINPGR')
!
!     *** user gives coord. corner point input grid               ***
!     *** Regular or curvilinear grid option  for version 30.21   ***
!
!     *** Type of grid -Regular(1) or Curvilinear(2)- is indicated***
!     *** by array IGTYPE(NUMGRD) in module SWCOMM2               ***
!
!     ------------------------------------------------------------------
!
!     INPUT    definition of input grids
!
! ========================================================================
!
!            | BOTtom   |
!            |          |
!            | WLEVel   |
!            |          |
!            | CURrent  |
!            |          |
!            | VX       |
!            |          |
!            | VY       |
! INPgrid  (<            >) &
!            | FRiction |
!            |          |
!            | WInd     |
!            |          |
!            | WX       |
!            |          |
!            | WY       |
!            |          |
!            | ASTD     |                                                 40.03
!            |          |
!            | NPLAnts  |                                                 40.55
!
!
!    | REGular [xpinp] [ypinp] [alpinp] [mxinp] [myinp] [dxinp] [dyinp] |
!    |                                                                  |
!   <  CURVilinear [stagrx] [stagry] [mxinp] [myinp]                     > &
!    |                                                                  |
!    | UNSTRUCtured                                                     |
!
!    (EXCeption  [excval])                                                 &
!
!                                        | -> SEC  |
!    (NONSTATionary [tbeginp] [deltinp] <     MIN   >  [tendinp])
!                                        |    HR   |
!                                        |    DAY  |
!
! =========================================================================
!
      CALL INKEYW ('STA', 'REG')                                          30.21
      IF (KEYWIS('CURV')) THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (4,                                                 32.02
     &       'Impossible option: 1D-simulation with curvilinear grid')    32.02
          RETURN
        ENDIF                                                             32.02
        IF (OPTG.EQ.5) THEN                                               40.80
           CALL MSGERR (4,                                                40.80
     &   'Curvilinear input grid cannot be used with unstructured grid')  40.80
           RETURN                                                         40.80
        ENDIF                                                             40.80
!
        IGTYPE(IGRID1) = 2                                                30.21
!       default values changed                                            30.60
        CALL INREAL ('STAGRX',STAGRX,'STA', 0.)                           30.60
        CALL INREAL ('STAGRY',STAGRY,'STA', 0.)                           30.60
        STAGX(IGRID1) = STAGRX                                            30.21
        STAGY(IGRID1) = STAGRY                                            30.21
        MXG(IGRID1) = MXG(IGRID1) - 1
        MYG(IGRID1) = MYG(IGRID1) - 1
!       default values changed                                            30.60
        CALL ININTG ('MXINP', MXG(IGRID1),'STA', MXC-1)                   30.60
        MXG(IGRID1) = MXG(IGRID1) + 1
        IF (ONED) THEN                                                    32.02
          CALL ININTG ('MYINP', MYG(IGRID1),'STA',0)                      30.70
          IF (MYG(IGRID1) .NE. 0) THEN                                    30.70
            CALL MSGERR (1, '1D-simulation: [myinp] set to zero !')       32.02
          ENDIF                                                           32.02
          MYG(IGRID1) = 0                                                 30.70
        ELSE
          CALL ININTG ('MYINP', MYG(IGRID1),'STA',MYC-1)                  30.80
        ENDIF                                                             32.02
        MYG(IGRID1) = MYG(IGRID1) + 1
      ELSEIF (KEYWIS('UNSTRUC')) THEN                                     40.80
!
        IF (ONED) THEN                                                    40.80
           CALL MSGERR (4,                                                40.80
     &            '1D-simulation cannot be done with unstructured grid')  40.80
           RETURN                                                         40.80
        ENDIF                                                             40.80
!
        IGTYPE(IGRID1) = 3                                                40.80
        MXG(IGRID1)    = nverts                                           40.80
        MYG(IGRID1)    = 1                                                40.80
        IF (IGRID2.GT.0) THEN
           IGTYPE(IGRID2) = 3                                             40.80
           MXG(IGRID2)    = nverts                                        40.80
           MYG(IGRID2)    = 1                                             40.80
        ENDIF
      ELSE
        CALL IGNORE('REG')
        IGTYPE(IGRID1) = 1                                                30.21
        IF (KSPHER .EQ. 0) THEN                                           40.13
          CALL READXY ('XPINP', 'YPINP', XPG(IGRID1), YPG(IGRID1),
     &                 'UNC', 0., 0.)                                     30.20
          CALL INREAL ('ALPINP',ALPG(IGRID1),'UNC',0.)                    30.20
        ELSE                                                              40.13
!         [xpinp] and [ypinp] are required in case of spherical coordinates
          CALL READXY ('XPINP', 'YPINP', XPG(IGRID1), YPG(IGRID1),
     &                 'REQ', 0., 0.)                                     40.13
          CALL INREAL ('ALPINP',ALPG(IGRID1),'STA',0.)                    40.13
        ENDIF                                                             40.13
!       ---------  ALPG(IGRID1) is always between -PI and PI   ---------
        ALTMP = ALPG(IGRID1)/360.
        ALPG(IGRID1)  = PI2 * (ALTMP - NINT(ALTMP))
        COSPG(IGRID1) = COS(ALPG(IGRID1))
        SINPG(IGRID1) = SIN(ALPG(IGRID1))
!       *** the user gives number of meshes in X resp. Y direction ***
!       *** the program uses the number of grid points             ***
        MXG(IGRID1) = MXG(IGRID1) - 1
        MYG(IGRID1) = MYG(IGRID1) - 1
        CALL ININTG ('MXINP', MXG(IGRID1),'RQI',-1)                       30.20
        IF (ONED) THEN                                                    32.02
          CALL ININTG ('MYINP', MYG(IGRID1),'STA',0)                      30.70
          IF (MYG(IGRID1) .NE. 0) THEN                                    32.02
            CALL MSGERR (1, '1D-simulation: [myinp] set to zero !')       32.02
            MYG(IGRID1) = 0                                               32.02
          ENDIF                                                           32.02
        ELSE
          CALL ININTG ('MYINP', MYG(IGRID1),'RQI',-1)                     30.20
        ENDIF                                                             32.02
        MXG(IGRID1) = MXG(IGRID1) + 1
        MYG(IGRID1) = MYG(IGRID1) + 1
!
        CALL INREAL ('DXINP',DXG(IGRID1),'RQI',0.)                        30.20
        CALL INREAL ('DYINP',DYG(IGRID1),'STA',DXG(IGRID1))               30.20
!
      ENDIF
!
!     exception values for input variables                                30.60
!
      CALL INKEYW ('STA', ' ')                                            30.60
      IF (KEYWIS ('EXC')) THEN                                            30.60
        CALL INREAL ('EXCVAL', EXCFLD(IGRID1), 'REQ', 0.)                 30.60
        IF (IGRID2.GT.0) EXCFLD(IGRID2) = EXCFLD(IGRID1)                  40.80
      ENDIF                                                               30.60
!
      LEDS(IGRID1) = 1
      IF (IGRID2.GT.0) LEDS(IGRID2) = 1                                   40.80
!     next lines to process the option NONSTAT for wind, currents         30.70
!     and waterlevel                                                VER.  30.00
!     SUBR. INCTIM reads a time string  and gives a time from given
!     reference day
      CALL INKEYW ('STA', ' ')
      IF (KEYWIS ('NONSTAT')) THEN
        IF (NSTATM.EQ.0) CALL MSGERR (3,
     &  'keyword NONSTAT not allowed in stationary mode')
        NSTATM = 1
        IF (IGRID1.EQ.1 .OR. IGRID1.EQ.8) CALL MSGERR (2,
     &        'nonstationary input field not allowed in this case')
        CALL INCTIM (ITMOPT,'TBEGINP',IFLBEG(IGRID1),'REQ',0.)            40.00
        CALL ININTV ('DELTINP', IFLINT(IGRID1), 'REQ', 0.)                40.00
        CALL INCTIM (ITMOPT,'TENDINP',IFLEND(IGRID1),'STA',1.E20)         40.00
        IFLDYN(IGRID1) = 1                                                40.00
        IFLTIM(IGRID1) = IFLBEG(IGRID1)                                   40.00
        IF (IGRID2 .GT. 0) THEN
          IFLBEG(IGRID2) = IFLBEG(IGRID1)                                 40.00
          IFLINT(IGRID2) = IFLINT(IGRID1)                                 40.00
          IFLEND(IGRID2) = IFLEND(IGRID1)                                 40.00
          IFLDYN(IGRID2) = IFLDYN(IGRID1)                                 40.00
          IFLTIM(IGRID2) = IFLTIM(IGRID1)                                 40.00
        ENDIF
        IF (IFLEND(IGRID1).LT.0.9E20) THEN
           IF ( MOD(IFLEND(IGRID1)-IFLBEG(IGRID1),IFLINT(IGRID1)) >
     &                                         0.01*IFLINT(IGRID1) .AND.
     &          MOD(IFLEND(IGRID1)-IFLBEG(IGRID1),IFLINT(IGRID1)) <
     &                                         0.99*IFLINT(IGRID1) )
     &        CALL MSGERR (1,
     &             '[deltinp] is not a fraction of the period')
        ENDIF
      ENDIF
      IF ((MXG(IGRID1).EQ.0).OR.(MYG(IGRID1).EQ.0)) GOTO 100
      IF (IGTYPE(IGRID1).EQ.3) GOTO 100                                   40.80
!
!     ***** input grid is included in output data                 *****
!     ***** reference point coincides with origin of output frame *****
      ALLOCATE(OPSTMP)                                                    40.31
      OPSTMP%PSNAME = SNAMEG                                              40.31
      ALLOCATE(OPSTMP%XP(0))                                              40.31
      ALLOCATE(OPSTMP%YP(0))                                              40.31
!
      IF (IGTYPE(IGRID1) .EQ. 1) THEN                                     30.21
        OPSTMP%PSTYPE = 'F'                                               40.31
        XQLEN          = (MXG(IGRID1)-1)*DXG(IGRID1)
        OPSTMP%OPR(3) = XQLEN                                             40.31
!
        IF (ONED) THEN                                                    32.02
          YQLEN          = XQLEN                                          32.02
        ELSE                                                              32.02
          YQLEN          = (MYG(IGRID1)-1)*DYG(IGRID1)
        ENDIF                                                             32.02
!
        OPSTMP%OPR(4) = YQLEN                                             40.31
        OPSTMP%OPR(1) = XPG(IGRID1)                                       40.31
        OPSTMP%OPR(2) = YPG(IGRID1)                                       40.31
        OPSTMP%OPR(5) = ALPG(IGRID1)                                      40.31
        OPSTMP%OPI(1) = MXG(IGRID1)                                       40.31
        OPSTMP%OPI(2) = MYG(IGRID1)                                       40.31
      ELSE
        OPSTMP%PSTYPE = 'H'                                               40.31
        OPSTMP%OPR(1) = FLOAT(MXG(IGRID1)-1)                              40.31
        OPSTMP%OPR(2) = FLOAT(MYG(IGRID1)-1)                              40.31
        OPSTMP%OPR(3) = 0.                                                40.31
        OPSTMP%OPR(4) = 0.                                                40.31
        OPSTMP%OPR(5) = 0.                                                40.31
        OPSTMP%OPI(1) = MXG(IGRID1)                                       40.31
        OPSTMP%OPI(2) = MYG(IGRID1)                                       40.31
!       *** Because the input grid is staggered should have    ***
!       *** same length in X and Y than the computational grid ***
        IF (MCGRD.GT.1) THEN
          XQLEN = XCGMAX - XCGMIN                                         40.00
          YQLEN = YCGMAX - YCGMIN
        ENDIF
      ENDIF
      IF (ITEST .GE. 50 .OR. INTES .GE. 5) THEN
        IF (OPSTMP%PSTYPE .EQ. 'F') THEN
          WRITE(PRINTF,6021)IGRID1,'F',XQLEN,YQLEN,
     &    XPG(IGRID1),YPG(IGRID1),ALPG(IGRID1),MXG(IGRID1),MYG(IGRID1)
        ELSE
          WRITE(PRINTF,6022)IGRID1,OPSTMP%PSTYPE,MXG(IGRID1)-1,
     &    MYG(IGRID1)-1,0,0,0.,MXG(IGRID1),MYG(IGRID1)
        ENDIF
 6021   FORMAT (' INP GRID PARAMETERS: ',/,
     &          'IGRID , FRAMTYPE,XLENFR ,YLENFR  XPFR  YPFR    ALPFR',
     &          '    MXFR MYFR',/,I2,1X,A,4X,
     &           2(1X,E8.3), 2(1X,E10.3), F7.3, 2(1X,I4))                 40.31 30.60
 6022   FORMAT (' INP GRID PARAMETERS: ',/,
     &          'IGRID, FRAMTYPE ,XMAXFR ,YMAXFR XMINFR YMINFR  ALPFR',
     &          '    MXFR MYFR',/,
     &          I3,8X,A,1X,4(3X,I4),5X,E8.3,2(1X,I4))
      ENDIF                                                               30.21
      NULLIFY(OPSTMP%NEXTOPS)                                             40.31
      IF ( .NOT.LOPS ) THEN                                               40.31
         FOPS = OPSTMP                                                    40.31
         COPS => FOPS                                                     40.31
         LOPS = .TRUE.                                                    40.31
      ELSE                                                                40.31
         COPS%NEXTOPS => OPSTMP                                           40.31
         COPS => OPSTMP                                                   40.31
      END IF                                                              40.31
!
      IF (IGRID2.GT.0) THEN                                               10.26
        XPG(IGRID2)   = XPG(IGRID1)
        YPG(IGRID2)   = YPG(IGRID1)
        ALPG(IGRID2)  = ALPG(IGRID1)
        COSPG(IGRID2) = COSPG(IGRID1)
        SINPG(IGRID2) = SINPG(IGRID1)
        DXG(IGRID2)   = DXG(IGRID1)
        DYG(IGRID2)   = DYG(IGRID1)
        MXG(IGRID2)   = MXG(IGRID1)
        MYG(IGRID2)   = MYG(IGRID1)
        STAGX(IGRID2) = STAGX(IGRID1)                                     40.80
        STAGY(IGRID2) = STAGY(IGRID1)                                     40.80
        IGTYPE(IGRID2)= IGTYPE(IGRID1)                                    30.51
      ENDIF                                                               10.26
!
      DO 80 IGRID = IGRID1+1, NUMGRD
        IF (LEDS(IGRID).EQ.0) THEN
          XPG(IGRID)   = XPG(IGRID1)
          YPG(IGRID)   = YPG(IGRID1)
          ALPG(IGRID)  = ALPG(IGRID1)
          COSPG(IGRID) = COSPG(IGRID1)
          SINPG(IGRID) = SINPG(IGRID1)
          DXG(IGRID)   = DXG(IGRID1)
          DYG(IGRID)   = DYG(IGRID1)
          MXG(IGRID)   = MXG(IGRID1)
          MYG(IGRID)   = MYG(IGRID1)
          IGTYPE(IGRID)= IGTYPE(IGRID1)                                   30.52
          LEDS(IGRID)  = 1
        ENDIF
  80  CONTINUE
!
 100  RETURN
!     end of subroutine SINPGR
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SREDEP ( LWINDR, LWINDM ,LOGCOM )                        40.31 40.02
!                                                                      *
!***********************************************************************
!
      USE TIMECOMM                                                        40.41
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_GENARR                                                        40.31
      USE SwanGriddata                                                    40.80
!ADC      USE Couple2Swan, ONLY: ADCIRC_ETA2 => SWAN_ETA2,                    41.20
!ADC     &                       ADCIRC_UU2 => SWAN_UU2,                      41.20
!ADC     &                       ADCIRC_VV2 => SWAN_VV2,                      41.20
!ADC     &                       ADCIRC_WX2 => SWAN_WX2,                      41.20
!ADC     &                       ADCIRC_WY2 => SWAN_WY2,                      41.20
!ADC     &                       COUPCUR, COUPWIND, COUPWLV                   41.20
!ADC     &                      ,ADCIRC_Z0 => SWAN_Z0,                        41.20
!ADC     &                       COUPFRIC                                     41.20
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
!     34.01: Jeroen Adema
!     40.02: IJsbrand Haagsma
!     40.04: Annette Kieftenburg
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.55: Marcel Zijlema
!     41.20: Casey Dietrich
!
!  1. Updates
!
!     20.05, Jan. 94: new pool
!     20.67, Dec. 95: call of REPARM is modified, VFAC is read in SREDEP itself
!     30.72, Nov. 97: Changed position of label 18 in block IF, as suggested by
!                     Richard Gorman
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.00, Jan. 98: calls of REPARM and INAR2D changed
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Oct. 00: Avoided real/int conflict by introducing RPOOL
!     40.31, Oct. 03: remove POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.55, Dec. 05: introducing vegetation model
!     41.20, Mar. 10: extension to tightly coupled ADCIRC+SWAN model
!
!  2. Purpose
!
!     Reading of depths and/or currents
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     LWINDR
!     LWINDM
!     LOGCOM
!
!  8. Subroutines used
!
!     INAR2D
!     INKEYW
!     KEYWIS
!     MSGERR
!     OTAR2D
!     REPARM (all Ocean Pack)
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWREAD
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     NDS   = unit reference number
!     DFORM = format string
!     VFAC  = multiplication factor for data to be read
!     DESCR = string used in heading of datafile
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     Call INKEYW to read keyword from user input
!     If bottom must be read (command BOT), then
!         If no current is read (LEDS < 2), then
!             LEDS = 1 (indicating INP BOT command has been given)        40.04
!             LEDS = 2 (indicating bottom is read)                        40.04
!         Else
!             LEDS = 3 (indicating bottom and current are read)
!         ------------------------------------------------------------
!         Call REPARM to process the rest of the command
!         If depthfile is standard Ocean Pack file (IDLA > 100), then
!             Call OTAR2D to read bottom levels into array IDEB
!         Else
!             Call INAR2D to read bottom levels inyo array IDEB
!         ------------------------------------------------------------
!     Elseif current must be read (command CUR), then
!         Switch for current is on (ICUR = 1)
!         If no current is read (LEDS < 2), then
!             LEDS = LEDS + 2 (indicating current is read)
!         ------------------------------------------------------------
!         Call REPARM to read filename and file organisation
!         If depthfile is standard Ocean Pack file (IDLA > 100), then
!             Call OTAR2D to read current components into arrays IUXB
!               and IUYB
!         Else
!             Call INAR2D to read current components into arrays IUXB
!               and IUYB
!         ------------------------------------------------------------
!     Else
!         Call MSGERR to generate an error message
!     ----------------------------------------------------------------
!
! 13. Source text
!
      REAL, ALLOCATABLE :: TARR(:)                                        40.31
      LOGICAL    KEYWIS, VECTOR, LOCAL                                    40.95 30.00
      LOGICAL    LOGCOM(6)                                                30.21
      SAVE  IENT
      DATA  IENT/0/
      CALL STRACE (IENT,'SREDEP')
!
!     ***** read instructions of the user *****
!
!   ============================================================================
!
!   READinp    BOTtom/WLevel/CURrent/FRiction/WInd/COORdinates               &
!              NPLAnts                                                       &
!        [fac]  / 'fname1'        \
!               \ SERIES 'fname2' /  [idla] [nhedf] ([nhedt]) (nhedvec])     &
!        FREE / FORMAT 'form' / [idfm] / UNFORMATTED
!
!   ============================================================================
!
      CALL INKEYW ('REQ',' ')
      IGR2 = 0
      IF (KEYWIS ('BOT')) THEN
        IGR1      = 1
        LOGCOM(3) = .TRUE.                                                30.21
      ELSE IF (KEYWIS ('CUR')) THEN
        IGR1 = 2
        IGR2 = 3
        ICUR = 1
!ADC      ELSE IF (KEYWIS ('ADCCUR')) THEN                                    41.20
!ADC!       added this possibility for coupling to ADCIRC currents
!ADC        IGR1 = 2
!ADC        IGR2 = 3
!ADC        ICUR = 1
!ADC        COUPCUR = .TRUE.
!ADC        IFLFAC(IGR1) = 1.0
!ADC        IFLNDF(IGR1) = 0
!ADC        IFLNDS(IGR1) = 22
!ADC        IFLIDL(IGR1) = 1
!ADC        IFLIFM(IGR1) = 2
!ADC        IFLFRM(IGR1) = '(12X,E15.10)'
!ADC        IFLNHF(IGR1) = 2
!ADC        IFLDYN(IGR1) = 1
!ADC        IFLNHD(IGR1) = 1
!ADC        VECTOR = .TRUE.
!ADC        NHEDC = 0
!ADC        IFLFAC(IGR2) = IFLFAC(IGR1)
!ADC        IFLNDF(IGR2) = IFLNDF(IGR1)
!ADC        IFLNDS(IGR2) = IFLNDS(IGR1)
!ADC        IFLIDL(IGR2) = IFLIDL(IGR1)
!ADC        IFLIFM(IGR2) = IFLIFM(IGR1)
!ADC        IFLFRM(IGR2) = IFLFRM(IGR1)
!ADC        IFLNHF(IGR2) = IFLNHF(IGR1)
!ADC        IFLNHD(IGR2) = NHEDC
!ADC        IF (.NOT.ALLOCATED(TARR)) THEN
!ADC           ALLOCATE( TARR(1:nverts) )
!ADC           TARR = 0.
!ADC        ENDIF
!ADC        DO I=1,nverts
!ADC           TARR(I) = REAL(ADCIRC_UU2(I,1))
!ADC        ENDDO
!ADC        GOTO 55
      ELSE IF (KEYWIS ('FR')) THEN
        IGR1   = 4
        VARFR  = .TRUE.
        MCMVAR = MCMVAR + 2                                               40.00
        JFRC2  = MCMVAR - 1                                               40.00
        JFRC3  = MCMVAR                                                   40.00
        ALOCMP = .TRUE.                                                   40.97
!ADC      ELSE IF (KEYWIS ('ADCFRIC')) THEN                                   41.20
!ADC!       added this possibility for coupling Madsen friction lengths
!ADC!       based on ADCIRC Manning's n values
!ADC        IGR1 = 4
!ADC        VARFR = .TRUE.
!ADC        MCMVAR = MCMVAR + 2
!ADC        JFRC2  = MCMVAR - 1
!ADC        JFRC3  = MCMVAR
!ADC        ALOCMP = .TRUE.
!ADC        COUPFRIC = .TRUE.
!ADC        IFLFAC(IGR1) = 1.0
!ADC        IFLNDF(IGR1) = 0
!ADC        IFLNDS(IGR1) = 23
!ADC        IFLIDL(IGR1) = 1
!ADC        IFLIFM(IGR1) = 2
!ADC        IFLFRM(IGR1) = '(12X,E20.10)'
!ADC        IFLNHF(IGR1) = 2
!ADC        IFLDYN(IGR1) = 1
!ADC        IFLNHD(IGR1) = 1
!ADC        VECTOR = .FALSE.
!ADC        NHEDC = 0
!ADC        IF (.NOT.ALLOCATED(TARR)) THEN
!ADC           ALLOCATE( TARR(1:nverts) )
!ADC           TARR = 0.
!ADC        ENDIF
!ADC        DO I=1,nverts
!ADC           TARR(I) = REAL(ADCIRC_Z0(I,1))
!ADC        ENDDO
!ADC        GOTO 55
      ELSE IF (KEYWIS ('WI')) THEN
        LWINDR = 2                                                        30.10
        IWIND  = LWINDM                                                   30.10
        IGR1   = 5
        IGR2   = 6
        VARWI  = .TRUE.
!ADC      ELSE IF (KEYWIS ('ADCWIND')) THEN                                   41.20
!ADC!       added this possibility for coupling to ADCIRC wind speeds
!ADC        LWINDR = 2
!ADC        IWIND  = LWINDM
!ADC        IGR1 = 5
!ADC        IGR2 = 6
!ADC        VARWI = .TRUE.
!ADC        COUPWIND = .TRUE.
!ADC        IFLFAC(IGR1) = 1.0
!ADC        IFLNDF(IGR1) = 0
!ADC        IFLNDS(IGR1) = 22
!ADC        IFLIDL(IGR1) = 1
!ADC        IFLIFM(IGR1) = 2
!ADC        IFLFRM(IGR1) = '(12X,E15.10)'
!ADC        IFLNHF(IGR1) = 2
!ADC        IFLDYN(IGR1) = 1
!ADC        IFLNHD(IGR1) = 1
!ADC        VECTOR = .TRUE.
!ADC        NHEDC = 0
!ADC        IFLFAC(IGR2) = IFLFAC(IGR1)
!ADC        IFLNDF(IGR2) = IFLNDF(IGR1)
!ADC        IFLNDS(IGR2) = IFLNDS(IGR1)
!ADC        IFLIDL(IGR2) = IFLIDL(IGR1)
!ADC        IFLIFM(IGR2) = IFLIFM(IGR1)
!ADC        IFLFRM(IGR2) = IFLFRM(IGR1)
!ADC        IFLNHF(IGR2) = IFLNHF(IGR1)
!ADC        IFLNHD(IGR2) = NHEDC
!ADC        IF (.NOT.ALLOCATED(TARR)) THEN
!ADC           ALLOCATE( TARR(1:nverts) )
!ADC           TARR = 0.
!ADC        ENDIF
!ADC        DO I=1,nverts
!ADC           TARR(I) = REAL(ADCIRC_WX2(I,1))
!ADC        ENDDO
!ADC        GOTO 55
      ELSE IF (KEYWIS ('WL')) THEN                                        20.38
        IGR1   = 7
        VARWLV = .TRUE.                                                   20.38
!ADC      ELSE IF (KEYWIS ('ADCWL')) THEN                                     41.20
!ADC!       added this possibility for coupling to ADCIRC water levels
!ADC        IGR1 = 7
!ADC        VARWLV = .TRUE.
!ADC        COUPWLV = .TRUE.
!ADC        IFLFAC(IGR1) = 1.0
!ADC        IFLNDF(IGR1) = 0
!ADC        IFLNDS(IGR1) = 23
!ADC        IFLIDL(IGR1) = 1
!ADC        IFLIFM(IGR1) = 2
!ADC        IFLFRM(IGR1) = '(12X,E20.10)'
!ADC        IFLNHF(IGR1) = 2
!ADC        IFLDYN(IGR1) = 1
!ADC        IFLNHD(IGR1) = 1
!ADC        VECTOR = .FALSE.
!ADC        NHEDC = 0
!ADC        IF (.NOT.ALLOCATED(TARR)) THEN
!ADC           ALLOCATE( TARR(1:nverts) )
!ADC           TARR = 0.
!ADC        ENDIF
!ADC        DO I=1,nverts
!ADC           TARR(I) = REAL(ADCIRC_ETA2(I,1))
!ADC        ENDDO
!ADC        GOTO 55
      ELSE IF (KEYWIS ('COOR')) THEN                                      30.21
        IGR1   = 8                                                        30.21
        IGR2   = 9                                                        30.21
        LOGCOM(4) = .TRUE.                                                30.21
!       *** Next lines because there is no information for READ COORD  ***
!       *** command INPGRID was not used to read coordinates           ***
        MXG(IGR1) = MXC
        MYG(IGR1) = MYC
        MXG(IGR2) = MXC
        MYG(IGR2) = MYC
      ELSE IF (KEYWIS ('ASTD')) THEN                                      40.03
!       air-sea temperature difference
        IGR1   = 10
        VARAST = .TRUE.                                                   40.03
        IF (JASTD2.LE.1) THEN                                             40.03
          MCMVAR = MCMVAR + 2                                             40.03
          JASTD2 = MCMVAR - 1                                             40.03
          JASTD3 = MCMVAR                                                 40.03
          ALOCMP = .TRUE.                                                 40.97
        ENDIF                                                             40.03
      ELSE IF (KEYWIS ('NPLA')) THEN                                      40.55
!       number of plants per square meter
        IGR1   = 11
        VARNPL = .TRUE.                                                   40.55
        IF (JNPLA2.LE.1) THEN                                             40.55
          MCMVAR = MCMVAR + 2                                             40.55
          JNPLA2 = MCMVAR - 1                                             40.55
          JNPLA3 = MCMVAR                                                 40.55
          ALOCMP = .TRUE.                                                 40.97
        ENDIF                                                             40.55
      ELSE
        CALL  WRNKEY
      ENDIF
!
!     read multiplication factor
!
      CALL INREAL ('FAC', IFLFAC(IGR1), 'STA', 1.)
!
      IF (IGR2.GT.0) THEN
        VECTOR = .TRUE.
        IFLFAC(IGR2) = IFLFAC(IGR1)
      ELSE
        VECTOR = .FALSE.
      ENDIF
!
      LOCAL = .FALSE.                                                     40.95
!PUN      IF ( IGTYPE(IGR1).EQ.3 ) THEN                                       40.95
!PUN         LOCAL = .TRUE.                                                   40.95
!PUN      ELSE                                                                40.95
!PUN         LOCAL = .FALSE.                                                  40.95
!PUN      ENDIF                                                               40.95
!
      CALL REPARM (IFLNDF(IGR1), IFLNDS(IGR1), IFLIDL(IGR1),              40.00
     &             IFLIFM(IGR1), IFLFRM(IGR1), IFLNHF(IGR1),              40.00
     &             IFLDYN(IGR1), IFLNHD(IGR1), VECTOR, LOCAL, NHEDC)      40.95 40.00
      IF (STPNOW()) RETURN                                                34.01
!
      IF (ITEST.GE.60) WRITE (PRTEST, 35) IGR1, IFLNDF(IGR1),
     &    IFLNDS(IGR1), IFLIDL(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),         40.00
     &    IFLNHF(IGR1), IFLDYN(IGR1), IFLNHD(IGR1), VECTOR, NHEDC
  35  FORMAT (' Reading parameters: ', 5I4, A, /, 12X,I5,I5,I5,L2,I5)     40.00
!
      IFLNHD(IGR1) = IFLNHD(IGR1) + NHEDC
      IF (IGR2.GT.0) THEN
        IFLNDF(IGR2) = IFLNDF(IGR1)
        IFLNDS(IGR2) = IFLNDS(IGR1)
        IFLIDL(IGR2) = IFLIDL(IGR1)
        IFLIFM(IGR2) = IFLIFM(IGR1)
        IFLFRM(IGR2) = IFLFRM(IGR1)
        IFLNHF(IGR2) = IFLNHF(IGR1)
        IFLNHD(IGR2) = NHEDC
      ENDIF
!
      IF (LEDS(IGR1).EQ.0 .AND. IGR1 .NE. 8) THEN                         30.21
         CALL MSGERR (2, 'Input grid not given')
         RETURN
      ENDIF

      IF ( IGR1.EQ.1 .AND. grid_generator.EQ.meth_adcirc ) THEN           40.80
!        bottom topography will be taken from fort.14
         CALL MSGERR(1,'depth will be taken from grid file fort.14 ')
         IGTYPE(1) = 3
         LEDS(1)   = 2
         RETURN
      ENDIF

      IF (.NOT.ALLOCATED(TARR)) THEN                                      40.41
         ALLOCATE( TARR(MXG(IGR1)*MYG(IGR1)) )                            40.41
         TARR = 0.                                                        40.41
      END IF                                                              40.41
      CALL INAR2D( TARR        , MXG(IGR1), MYG(IGR1), IFLNDF(IGR1),      40.31 40.02
     &             IFLNDS(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),
     &             IFLIDL(IGR1), IFLFAC(IGR1),
     &             IFLNHD(IGR1), IFLNHF(IGR1))
      IF (STPNOW()) RETURN                                                34.01
  55  CONTINUE                                                            41.20
      IF (IGR1.EQ.1) THEN                                                 40.31
         IF (.NOT.ALLOCATED(DEPTH)) ALLOCATE(DEPTH(MXG(IGR1)*MYG(IGR1)))  40.31
         CALL SWCOPR( TARR, DEPTH, MXG(IGR1)*MYG(IGR1) )                  40.31
      ELSE IF (IGR1.EQ.2) THEN                                            40.31
         IF (.NOT.ALLOCATED(UXB)) ALLOCATE(UXB(MXG(IGR1)*MYG(IGR1)))      40.31
         CALL SWCOPR( TARR, UXB, MXG(IGR1)*MYG(IGR1) )                    40.31
      ELSE IF (IGR1.EQ.4) THEN                                            40.31
         IF (.NOT.ALLOCATED(FRIC)) ALLOCATE(FRIC(MXG(IGR1)*MYG(IGR1)))    40.31
         CALL SWCOPR( TARR, FRIC, MXG(IGR1)*MYG(IGR1) )                   40.31
      ELSE IF (IGR1.EQ.5) THEN                                            40.31
         IF (.NOT.ALLOCATED(WXI)) ALLOCATE(WXI(MXG(IGR1)*MYG(IGR1)))      40.31
         CALL SWCOPR( TARR, WXI, MXG(IGR1)*MYG(IGR1) )                    40.31
      ELSE IF (IGR1.EQ.7) THEN                                            40.31
         IF (.NOT.ALLOCATED(WLEVL)) ALLOCATE(WLEVL(MXG(IGR1)*MYG(IGR1)))  40.31
         CALL SWCOPR( TARR, WLEVL, MXG(IGR1)*MYG(IGR1) )                  40.31
      ELSE IF (IGR1.EQ.8) THEN                                            40.31
         CALL SWCOPR( TARR, XCGRID, MXG(IGR1)*MYG(IGR1) )                 40.31
      ELSE IF (IGR1.EQ.10) THEN                                           40.31
         IF (.NOT.ALLOCATED(ASTDF)) ALLOCATE(ASTDF(MXG(IGR1)*MYG(IGR1)))  40.31
         CALL SWCOPR( TARR, ASTDF, MXG(IGR1)*MYG(IGR1) )                  40.31
      ELSE IF (IGR1.EQ.11) THEN                                           40.55
         IF (.NOT.ALLOCATED(NPLAF)) ALLOCATE(NPLAF(MXG(IGR1)*MYG(IGR1)))  40.55
         CALL SWCOPR( TARR, NPLAF, MXG(IGR1)*MYG(IGR1) )                  40.55
      END IF                                                              40.31
      DEALLOCATE(TARR)                                                    40.31
      IF (IGR2.GT.0) THEN

        IF (.NOT.ALLOCATED(TARR)) THEN                                    40.41
           ALLOCATE( TARR(MXG(IGR2)*MYG(IGR2)) )                          40.41
           TARR = 0.                                                      40.41
        END IF                                                            40.41
!ADC        IF ( IGR2.EQ.3 .AND. COUPCUR ) THEN                               41.20
!ADC           DO I = 1, nverts                                               41.20
!ADC              TARR(I) = REAL(ADCIRC_VV2(I,1))                             41.20
!ADC           ENDDO                                                          41.20
!ADC        ELSE IF ( IGR2.EQ.6 .AND. COUPWIND ) THEN                         41.20
!ADC           DO I = 1, nverts                                               41.20
!ADC              TARR(I) = REAL(ADCIRC_WY2(I,1))                             41.20
!ADC           ENDDO                                                          41.20
!ADC        ELSE                                                              41.20
        CALL INAR2D( TARR        , MXG(IGR2), MYG(IGR2), IFLNDF(IGR2),    40.31 40.02
     &               IFLNDS(IGR2), IFLIFM(IGR2), IFLFRM(IGR2),
     &               IFLIDL(IGR2), IFLFAC(IGR2),
     &               IFLNHD(IGR2), 0)
        IF (STPNOW()) RETURN                                              34.01
!ADC        ENDIF                                                             41.20
        IF (IGR2.EQ.3) THEN                                               40.31
           IF (.NOT.ALLOCATED(UYB)) ALLOCATE(UYB(MXG(IGR2)*MYG(IGR2)))    40.31
           CALL SWCOPR( TARR, UYB, MXG(IGR2)*MYG(IGR2) )                  40.31
        ELSE IF (IGR2.EQ.6) THEN                                          40.31
           IF (.NOT.ALLOCATED(WYI)) ALLOCATE(WYI(MXG(IGR2)*MYG(IGR2)))    40.31
           CALL SWCOPR( TARR, WYI, MXG(IGR2)*MYG(IGR2) )                  40.31
        ELSE IF (IGR2.EQ.9) THEN                                          40.31
           CALL SWCOPR( TARR, YCGRID, MXG(IGR2)*MYG(IGR2) )               40.31
        END IF                                                            40.31
        DEALLOCATE(TARR)                                                  40.31
      ENDIF
!     set time of reading
      IF (IFLDYN(IGR1) .EQ. 1) THEN
        IF (NSTATM.EQ.0) CALL MSGERR (2,                                  40.00
     &        'nonstationary input field requires MODE NONSTAT')          40.00
        IFLTIM(IGR1) = IFLBEG(IGR1)
        IF (IGR2.GT.0) IFLTIM(IGR2) = IFLTIM(IGR1)
      ENDIF
!
      LEDS(IGR1) = 2
      IF (VECTOR) LEDS(IGR2) = 2
!
      RETURN
! * end of subroutine SREDEP *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SSFILL (SPCSIG, SPCDIR)                                  30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_WCAP                                                          40.02
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
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     20.43         : Calculation of spectral directions added
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of SPCDIR
!     40.02, Sep. 00: Calculate powers of Sigma for later use
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Discretisation in frequency (sigma) and direction (theta)
!
!  3. Method
!
!     Create a logarithmic distribution in frequency between lowest and highest
!     frequencies and store them in the array SPCSIG.
!
!     Create a distribution in direction that does not coincide with the orientation
!     of the computational grid, calculate some geometric derivatives and store
!     it in the array SPCDIR
!
!  4. Argument variables
!
!   o SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
!   o SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!
!  5. SUBROUTINES CALLING
!
!       none
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ---
!
! 10. SOURCE TEXT
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SSFILL')
!
!     Allocate arrays in module M_WCAP                                    40.02
!
      IF (.NOT.ALLOCATED(SIGPOW)) ALLOCATE (SIGPOW(MSC,6))                40.02
!
!
!     distribution of spectral frequencies
!
!     FRINTF is the frequency integration factor (=df/f)                  20.35
      FRINTF = ALOG(SHIG/SLOW) / FLOAT(MSC-1)                             20.35
      SFAC   = EXP(FRINTF)                                                20.35
      FRINTH = SQRT(SFAC)
!     determine spectral frequencies (logarithmic distribution)
      SPCSIG(1) = SLOW                                                    30.72
      DO 10 IS = 2, MSC
         SPCSIG(IS) = SPCSIG(IS-1) * SFAC                                 30.72
  10  CONTINUE
!
!     Calculate powers of sigma and store in global array
!
      SIGPOW(:,1) = SPCSIG                                                40.02
      SIGPOW(:,2) = SPCSIG**2                                             40.02
      SIGPOW(:,3) = SPCSIG * SIGPOW(:,2)                                  40.02
      SIGPOW(:,4) = SPCSIG * SIGPOW(:,3)                                  40.02
      SIGPOW(:,5) = SPCSIG * SIGPOW(:,4)                                  40.02
      SIGPOW(:,6) = SPCSIG * SIGPOW(:,5)                                  40.02
!
!     distribution of spectral directions
!
      DO 20 ID = 1, MDC
         SPCDIR(ID,1) = SPDIR1 + FLOAT(ID-1)*DDIR                         20.43
!        if a direction coincides with a direction of the (regular)       40.00
!        computational grid it is slightly changed
         IF (OPTG.EQ.1) THEN                                              40.00
           IF (ABS(MODULO(SPCDIR(ID,1)-ALPC+0.25*PI ,0.5*PI)-0.25*PI)     40.03
     &         .LT. 1.E-6) THEN
             OLDDIR = SPCDIR(ID,1)
             SPCDIR(ID,1) = OLDDIR + 2.E-6                                40.03
             IF (ITEST.GE.50) WRITE (PRINTF, 24) ID,
     &                      OLDDIR*180./PI, SPCDIR(ID,1)*180./PI          40.03
  24         FORMAT (' Modified spectral direction', I4, 2(2X,F10.5))     40.03
           ENDIF
         ENDIF                                                            40.00
         SPCDIR(ID,2) = COS(SPCDIR(ID,1))                                 20.43
         SPCDIR(ID,3) = SIN(SPCDIR(ID,1))                                 20.43
         SPCDIR(ID,4) = SPCDIR(ID,2) **2                                  20.43
         SPCDIR(ID,5) = SPCDIR(ID,2) * SPCDIR(ID,3)                       20.43
         SPCDIR(ID,6) = SPCDIR(ID,3) **2                                  20.43
  20  CONTINUE
!
      RETURN
! * end of function SSFILL *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CGINIT (LOGCOM)                                          40.31 30.90
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
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
!     30.81: Annette Kieftenburg
!     30.90: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.04: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     40.00, June 98: new subroutine replacing code inside SWREAD
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.81, Nov. 98: Adjustment for 1-D case of new boundary conditions
!     34.01, Feb. 99: Introducing STPNOW
!     40.04, Aug. 00: adjusted argumentlist of CGBOUN
!     40.02, Oct. 00: Avoided real/int conflict by introducing replacing
!                     RPOOL for POOL in DPPUTR
!     40.30, Apr. 03: introduction distributed-memory approach using MPI
!     40.31, Dec. 03: removing POOL-mechanism and reconsidering this
!                     subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Jun. 07: extension to unstructured grids
!
!  2. PURPOSE
!
!     Initialise arrays for description of computational grid
!
!  3. METHOD
!
!
!  4. Argument variables
!
!     LOGCOM:
!
      LOGICAL LOGCOM(6)
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR                  40.41
!     IARR  :     help array                                              40.31
!     IENT  :     number of entries                                       40.31
!     INDX  :     index counter for global grid                           40.31
!     IX    :     loop counter                                            40.31
!     IY    :     loop counter                                            40.31
!     MCGRDL:     number of wet grid points in own subdomain              40.31
!     MSGSTR:     string to pass message to call MSGERR                   40.41
!
      INTEGER IENT, INDX, IX, IY, MCGRDL                                  40.31
      INTEGER ISTAT, IF1, IL1                                             40.41
      INTEGER, ALLOCATABLE :: IARR(:)                                     40.31
      CHARACTER*20 NUMSTR, CHARS(1)                                       40.41
      CHARACTER*80 MSGSTR                                                 40.41
!
!  8. SUBROUTINES CALLING
!
!     SWREAD
!
!  9. SUBROUTINES USED
!
!     CGBOUN           Determines boundary of computational region        40.31
!     MSGERR : Handles error messages according to severity               40.41
!     NUMSTR : Converts integer/real to string                            40.41
!     STRACE           Tracing routine for debugging                      40.31
!     SWDECOMP                                                            40.31
!     SWCOPI                                                              40.31
!TIMG!     SWTSTA                                                              40.31
!TIMG!     SWTSTO                                                              40.31
!     TXPBLA : Removes leading and trailing blanks in string              40.41
!
      LOGICAL STPNOW                                                      34.01
!
! 10. ERROR MESSAGES
!
!     ---
!
! 11. REMARKS
!
!     ---
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'CGINIT')
!
!     --- Calculation of computational grid dimension ***
!
      IF(.NOT.ALLOCATED(KGRPGL)) ALLOCATE(KGRPGL(MXC,MYC))                40.31
!
      CALL SWDIM ( KGRPGL, DEPTH, XCGRID, YCGRID )                        40.31 30.60
!
!     --- call CGBOUN to determine computational grid outline
!
      IF (ONED) THEN                                                      30.81
        IF(.NOT.ALLOCATED(KGRBGL)) ALLOCATE(KGRBGL(4))                    40.31
        CALL CGBOUN ( KGRPGL, KGRBGL )                                    40.31 40.04
      ELSE                                                                30.81
        IF(.NOT.ALLOCATED(IARR)) ALLOCATE(IARR(2*MCGRD))                  40.31
        CALL CGBOUN ( KGRPGL, IARR )                                      40.31 40.04
        IF ( NGRBND.GT.0 ) THEN                                           40.31
           IF(.NOT.ALLOCATED(KGRBGL)) ALLOCATE(KGRBGL(2*NGRBND))          40.31
           CALL SWCOPI(IARR,KGRBGL,2*NGRBND)                              40.31
        ELSE                                                              40.31
           IF(.NOT.ALLOCATED(KGRBGL)) ALLOCATE(KGRBGL(0))                 40.31
        END IF                                                            40.31
      END IF                                                              30.81
      NGRBGL = NGRBND                                                     40.30
!
!     --- Carry out domain decomposition meant for                        40.31
!         distributed-memory approach                                     40.31

!TIMG      CALL SWTSTA(211)                                                    40.31
      CALL SWDECOMP                                                       40.31
!TIMG      CALL SWTSTO(211)                                                    40.31
      IF (STPNOW()) RETURN                                                40.31

!TIMG      CALL SWTSTA(212)                                                    40.31

!     --- Create copy of parts of KGRPGL for each subdomain -> KGRPNT     40.31

      IF(.NOT.ALLOCATED(KGRPNT)) ALLOCATE(KGRPNT(MXC,MYC))                40.31
      KGRPNT = 1                                                          40.31
      MCGRDL = 1                                                          40.31
      DO IX = MXF, MXL                                                    40.31
         DO IY = MYF, MYL                                                 40.31
            INDX = KGRPGL(IX,IY)                                          40.31
            IF ( INDX.NE.1 ) THEN                                         40.31
               MCGRDL = MCGRDL + 1                                        40.31
               KGRPNT(IX-MXF+1,IY-MYF+1) = MCGRDL                         40.31
            END IF                                                        40.31
         END DO                                                           40.31
      END DO                                                              40.31

!     --- Create copies of parts of XCGRID, YCGRID for each subdomain     40.31

      IF (.NOT.ALLOCATED(XGRDGL)) ALLOCATE(XGRDGL(MXCGL,MYCGL))           40.31
      IF (.NOT.ALLOCATED(YGRDGL)) ALLOCATE(YGRDGL(MXCGL,MYCGL))           40.31
      XGRDGL = XCGRID                                                     40.31
      YGRDGL = YCGRID                                                     40.31
      DEALLOCATE(XCGRID,YCGRID)                                           40.31
      ALLOCATE(XCGRID(MXC,MYC))                                           40.31
      ALLOCATE(YCGRID(MXC,MYC))                                           40.31
      DO IX = MXF, MXL                                                    40.31
         DO IY = MYF, MYL                                                 40.31
            XCGRID(IX-MXF+1,IY-MYF+1) = XGRDGL(IX,IY)                     40.31
            YCGRID(IX-MXF+1,IY-MYF+1) = YGRDGL(IX,IY)                     40.31
         END DO                                                           40.31
      END DO                                                              40.31

!     --- Computation of XCLMIN, XCLMAX, YCLMIN, YCLMAX                   40.41

      XCLMIN =  1.E09
      YCLMIN =  1.E09
      XCLMAX = -1.E09
      YCLMAX = -1.E09
      DO IX = 1, MXC
         DO IY = 1, MYC
            IF (KGRPNT(IX,IY).GT.1) THEN
               IF (XCGRID(IX,IY).LT.XCLMIN) XCLMIN = XCGRID(IX,IY)
               IF (YCGRID(IX,IY).LT.YCLMIN) YCLMIN = YCGRID(IX,IY)
               IF (XCGRID(IX,IY).GT.XCLMAX) XCLMAX = XCGRID(IX,IY)
               IF (YCGRID(IX,IY).GT.YCLMAX) YCLMAX = YCGRID(IX,IY)
            END IF
         END DO
      END DO

!     --- Create copy of parts of KGRBGL for each subdomain -> KGRBND     40.31

      IF ( NGRBND.GT.0 ) THEN                                             40.31
         IF (PARLL) THEN
            IARR = 0                                                      40.31
            CALL CGBOUN ( KGRPNT, IARR )                                  40.31
            IF(.NOT.ALLOCATED(KGRBND)) ALLOCATE(KGRBND(2*NGRBND))         40.31
            IF (NGRBND.GT.0) CALL SWCOPI(IARR,KGRBND,2*NGRBND)            40.31
         ELSE                                                             40.31
            IF(.NOT.ALLOCATED(KGRBND)) ALLOCATE(KGRBND(2*NGRBND))         40.31
            KGRBND = KGRBGL                                               40.31
         END IF                                                           40.31
      ELSE                                                                40.31
         IF(.NOT.ALLOCATED(KGRBND)) ALLOCATE(KGRBND(0))                   40.31
      END IF                                                              40.31
      IF(ALLOCATED(IARR)) DEALLOCATE(IARR)                                40.31

!TIMG      CALL SWTSTO(212)                                                    40.31

      IF(.NOT.ALLOCATED(AC2)) ALLOCATE(AC2(MDC,MSC,MCGRD),STAT=ISTAT)     40.41 40.31
      IF ( ISTAT.NE.0 ) THEN                                              40.41
         CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                             40.41
         CALL TXPBLA(CHARS(1),IF1,IL1)                                    40.41
         MSGSTR = 'Allocation problem: array AC2 and return code is '//   40.41
     &            CHARS(1)(IF1:IL1)                                       40.41
         CALL MSGERR ( 4, MSGSTR )                                        40.41
         RETURN                                                           40.41
      END IF                                                              40.41
      AC2 = 0.                                                            40.31
      LOGCOM(6) = .TRUE.
!
!     Note: piece of code w.r.t. defining COMPGRID has been               40.31
!           moved to command CGRID in routine SWREAD!                     40.31
!
!     --- the following arrays for unstructured grids are allocated       40.80
!         as empty ones                                                   40.80
!
      IF ( .NOT.ALLOCATED(xcugrd) ) ALLOCATE(xcugrd(0))                   40.80
      IF ( .NOT.ALLOCATED(ycugrd) ) ALLOCATE(ycugrd(0))                   40.80
      IF ( .NOT.ALLOCATED( vmark) ) ALLOCATE( vmark(0))                   40.80
!
      RETURN
!     end of subroutine CGINIT
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWDIM ( KGRPNT, DEPTH, XCGRID, YCGRID )                  40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
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
!     30.60: Nico Booij
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.60, July 97: exception values for coordinates are introduced
!                     subroutine restructured
!     30.60, Aug. 97: correction KGRPNT
!     30.60, Aug. 97: test point introduced
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.03, Dec. 99: computation of XCGMIN, XCGMAX, YCGMIN, YCGMAX is now done
!                     for curvilinear and regular grids.
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     ---
!
!  3. Method
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!  5. SUBROUTINES CALLING
!
!       SWREAD
!
!  6. SUBROUTINES USED
!
!       SVALQI (SWAN/SER)
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       For all potential grid points do
!           Make grid address =1 (means invalid grid address)
!           Determine whether this point is a test point                  30.60
!           If coordinates are valid
!           Then If grid offset is not yet defined (LXOFFS = False)
!                Then make grid offset equal to coordinates of this point
!                     Make LXOFFS = True
!                Else Subtract offset from coordinates
!                --------------------------------------------------------
!                If bottom level is not an exception value
!                Then assign a valid grid address to this grid point
!                     increase MCGRD by 1
!       -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      INTEGER     KGRPNT(MXC,MYC)
      REAL        DEPTH(*)                                                40.00
      LOGICAL     EQREAL                                                  40.00
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SWDIM')
!
      DO 62 IX = 1, MXC                                                   30.72
        DO 61 IY = 1, MYC
!         at start make each grid point invalid
          KGRPNT(IX,IY) = 1
!
!         If coordinates are valid (excluding exception values)           30.60
!         Then If grid offset is not yet defined (LXOFFS = False)
!              Then make grid offset equal to coordinates of this point
!                   Make LXOFFS = True
!              Else Subtract offset from coordinates
!
          IF (EQREAL(XCGRID(IX,IY), EXCFLD(8))) THEN                      40.00
            IF (.NOT. EQREAL(YCGRID(IX,IY), EXCFLD(9))) THEN              40.00
              CALL MSGERR (2, 'incorrect grid coordinates')
              WRITE (PRINTF, 811) XCGRID(IX,IY), YCGRID(IX,IY)            30.72
 811          FORMAT (' X= ', E12.4, '  Y= ', E12.4)
            ENDIF
          ELSE
            IF (EQREAL(YCGRID(IX,IY), EXCFLD(9))) THEN                    40.00
              CALL MSGERR (2, 'incorrect grid coordinates')
              WRITE (PRINTF, 811) XCGRID(IX,IY), YCGRID(IX,IY)            30.72
            ELSE
              IF (OPTG.EQ.3) THEN                                         30.60
                IF (.NOT. LXOFFS) THEN
                  XOFFS  = XCGRID(IX,IY)                                  30.72
                  YOFFS  = YCGRID(IX,IY)                                  30.72
                  LXOFFS = .TRUE.
                  XCGRID(IX,IY) = 0.
                  YCGRID(IX,IY) = 0.
                ELSE
                  XCGRID(IX,IY) = REAL(XCGRID(IX,IY) - DBLE(XOFFS))       30.72
                  YCGRID(IX,IY) = REAL(YCGRID(IX,IY) - DBLE(YOFFS))       30.72
                ENDIF
              ENDIF
              XP = XCGRID(IX,IY)                                          30.72
              YP = YCGRID(IX,IY)                                          30.72
!               ***** compute bottom level *****
!
              DEP = SVALQI (XP, YP, 1, DEPTH, 1 ,IX ,IY)                  40.00
!
!             If bottom level is not an exception value
!             Then assign a valid grid address to this grid point
!                  increase MCGRD by 1
              IF (.NOT. EQREAL(DEP, EXCFLD(1))) THEN                      40.00
                MCGRD = MCGRD + 1                                         30.60
                KGRPNT(IX,IY) = MCGRD                                     30.60
              ENDIF
              IF (ITEST .GE. 250 .OR. INTES .GE. 30)                      40.31 30.60
     &          WRITE (PRINTF,30) IX, IY, XP, YP,
     &                    DEP, KGRPNT(IX,IY)                              30.60
  30          FORMAT(2(I3,1X),1X,3(F10.1,1X),5X,I5)
            ENDIF
          ENDIF                                                           30.60
  61    CONTINUE                                                          30.72
  62  CONTINUE                                                            30.72
!
      EXCFLD(8) = REAL(EXCFLD(8) - DBLE(XOFFS))
      EXCFLD(9) = REAL(EXCFLD(9) - DBLE(YOFFS))
!
      IF (MCGRD.LE.1) CALL MSGERR (3, 'No valid grid points found')       30.60
      IF (ITEST.GE.60) WRITE(PRINTF,*)
     &    ' Offset values in SWDIM:', XOFFS, YOFFS,                       30.60
     &    ' ; ', MCGRD-1, ' grid points'                                  30.60
!
!     check geometric validity of the grid (all meshes must have same
!     orientation when going around the mesh)
      CALL CVCHEK (KGRPNT, XCGRID, YCGRID)                                30.72
!
!     *** Computation of XCGMIN, XCGMAX, YCGMIN, YCGMAX ***               40.03
      XCGMIN =  1.E09                                                     40.00
      YCGMIN =  1.E09
      XCGMAX = -1.E09
      YCGMAX = -1.E09
      DO 60 IX = 1, MXC
        DO 59 IY = 1, MYC                                                 30.72
          IF (KGRPNT(IX,IY) .GT. 1) THEN
            IF (XCGRID(IX,IY) .LT. XCGMIN) XCGMIN = XCGRID(IX,IY)         40.00
            IF (YCGRID(IX,IY) .LT. YCGMIN) YCGMIN = YCGRID(IX,IY)
            IF (XCGRID(IX,IY) .GT. XCGMAX) XCGMAX = XCGRID(IX,IY)
            IF (YCGRID(IX,IY) .GT. YCGMAX) YCGMAX = YCGRID(IX,IY)
          ENDIF
 59     CONTINUE                                                          30.72
 60   CONTINUE                                                            30.72
!     *** Computation of xclen and yclen in a curvilinear grid case ***
      IF (OPTG .EQ. 3) THEN                                               40.80 40.03
        XCLEN = XCGMAX - XCGMIN
        YCLEN = YCGMAX - YCGMIN
      ENDIF
      IF (ITEST .GE. 100) THEN
        WRITE (PRINTF,*)' Min and Max X and Y from subr SWDIM',
     &         XCGMIN+XOFFS, XCGMAX+XOFFS, YCGMIN+YOFFS, YCGMAX+YOFFS     40.00
        WRITE (PRINTF,*)' Size in X and Y from subr SWDIM',
     &         XCLEN, YCLEN
      ENDIF
!
      RETURN
! * end of subroutine SWDIM *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CGBOUN (KGRPNT, KGRBND)                                  40.04
!                                                                      *
!***********************************************************************

      USE SWCOMM3                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31

      IMPLICIT NONE
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
!     40.00, 40.03, 40.13  Nico Booij
!     30.81  Annette Kieftenburg
!     40.04  Annette Kieftenburg
!     40.30  Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     New function for curvilinear version (ver. 40.00). May '98
!     30.81, Nov. 98: Adjustment for 1-D case of new boundary conditions
!     40.03, Dec. 99: 2d procedure modified; boundary can now consist also
!                     of diagonals
!     40.04, Aug. 00: prevented that boundary is not closed: several checks
!                     added
!                     argument list adjusted
!     40.13, July 01: 1-D procedure corrected
!     40.30, Apr. 01: introduction distributed-memory approach using MPI
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Determine array containing all points of (a) (closed) boundary/boundaries
!     within the computational grid
!
!  3. Method (updated...)
!
!     1D: go through all gridpoints (ascending and descending) to check
!         for validity
!         save first and last point in boundary array
!     2D: For all grid points:
!         0) find first point which is possibly on the boundary
!         1) check whether neighbour is a valid point and valid boundary point
!            save information that point is scanned
!         2) If so store point in boundary array and repeat from 1
!            Else go to next neighbour (if not scanned already
!                                       Else remove isolated point from
!                                            computational grid and start from 0)
!
!  4. Argument variables
!
!     KGRPNT    in-& output   indirect addresses for grid points
!     KGRBND    output        array containing all boundary points
!                             (+ 2 extra zeros as area separator
!                             for all separated areas)
!
      INTEGER KGRPNT(MXC,MYC), KGRBND(*)
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ICGRD         Counter for computational gridcells
!     IDIR          Direction number 1 = to the right
!                                    2 = upwards
!                                    3 = to the left
!                                    4 = downwards
!     IENT          Number of entries of this subroutine
!     IXC, IYC      X- and Y-index of point under consideration
!     IXNEW, IYNEW  X- and Y-index of point under consideration
!     IXOLD, IYOLD  X- and Y-index of point under consideration
!     KGRPNTNEW     Indirect addresses for grid points
!     MCGRDNEW      Number of points in computational grid after elimination
!                   of isolated and points that are part of 1D configurations
!                   or connections
!     MCGRDOLD      Number of points in computational grid before elimination
!                   of isolated and points that are part of 1D configurations
!                   or connections
!     WNP           Number of Wet Neighbouring Points
!     SCANND        Array containing information whether point is scanned
!
!
      INTEGER  ICGRD, IDIR, IENT, IXC, IXNEW, IXOLD
      INTEGER  IYC, IYNEW, IYOLD
      INTEGER  MCGRDNEW, MCGRDOLD, WNP
      INTEGER, ALLOCATABLE :: KGRPNTNEW(:,:), SCANND(:)
!
!  8. Subroutines used
!
!     Function PVALID
!     Function VALIDBP
!     STRACE
!     SEPARAREA
!     MSGERR
!
      LOGICAL   PVALID, VALIDBP
!
!  9. Subroutines calling
!
!     CGINIT
!
! 10. Error messages
!
! 11. Remarks
!
!     In case of 2D:
!     A boundary is scanned only once: 1D areas or connections between
!     areas are excluded
!
! 12. Structure
!
!     -----------------------------------------------------------------
!     If 1D: For all gridpoints (ascending)
!            If gridpoint is valid store it's value in KGRBND(1)
!            Else give error message
!            For all gridpoints (descending)
!            If gridpoint is valid store it's value in KGRBND(3)
!            Else give error message
!     Number of boundary outline points is 2
!     -----------------------------------------------------------------
!     Else:
!     Make number of boundary outline points = 0
!     For all computational grid points do
!         make SCANNED(ix,iy) = False
!     Save MCGRD
!     -----------------------------------------------------------------
!     For all computational grid points do
!         If (ix,iy) is a valid grid point
!         Then If point (ix,iy) is no valid boundary point or
!                 point (ix,iy) has neighbouring points that are all
!                               invalid boundary points
!              Then remove this point from computational grid
!                   SCANND(ix,iy) = True
!         If (ix,iy) is not scanned
!              If point (ix-1,iy) is not a valid grid point and
!                 (ix,iy) is a valid boundary point
!              Then Make ixold=ix; iyold=iy
!                   Increase number of boundary outline points by 1
!                   Store (ixold,iyold)
!                   Make SCANNED(ixold,iyold) = True
!                   Make idir=4
!                   If number of Wet Neighbouring points = 3
!                   possibly two areas should be separated
!                   Repeat
!                       Case idir=
!                       =1: Make ixnew=ixold+1; iynew=iyold
!                       =2: Make ixnew=ixold  ; iynew=iyold+1
!                       =3: Make ixnew=ixold-1; iynew=iyold
!                       =4: Make ixnew=ixold  ; iynew=iyold-1
!                       -----------------------------------------------
!                       If (ixnew,iynew) is a valid grid point and
!                          valid boundary point
!                       Then If SCANNED(ixnew,iynew)
!                            Then exit from repeat
!                            Else
!                              Increase number of boundary outline points by 1
!                              Store (ixnew,iynew) in KGRBND array
!                              Make SCANNED(ixnew,iynew) = True
!                              Make ixold=ixnew; iyold=iynew
!                              Make idir = idir-1
!                              If number of Wet Neighbouring points = 3
!                              possibly two areas should be separated
!                              If idir=0
!                              Then Make idir=4
!                        Else Make idir = idir+1
!                          If (ixnew,iynew) is an invalid boundary point
!                             (ixnew,iynew) is a valid grid point
!                             remove point from computational grid
!                             SCANND(ixnew,iynew) = True
!                          If idir=5
!                          Then Make idir=1
!                   ---------------------------------------------------
!                   Increase number of boundary outline points by 1
!                   Store (0,0)      {separation between outlines}
!     -----------------------------------------------------------------
!     If MCGRD has changed
!     Then count number of valid gridpoints
!          store  indirect addressing number in new array
!       give MCGRD new value
!       write new information to old array
!     -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT,'CGBOUN')
!
      IF (ONED) THEN                                                      30.81
        KGRBND(1) = -1                                                    40.13 30.81
        KGRBND(2) = 1                                                     40.13 30.81
        DO IXC = 1, MXC                                                   30.81
          IF (KGRPNT(IXC,1).GT.1) THEN                                    30.81
            KGRBND(1) = IXC                                               30.81
            GOTO 81                                                       30.81
          END IF                                                          30.81
  81    ENDDO                                                             30.81
        IF (KGRBND(1) .LT. 0) THEN                                        40.13 30.81
          CALL MSGERR(3,'No valid gridpoint defined')                     40.13 30.81
        END IF                                                            40.13 30.81
        KGRBND(3) = -1                                                    40.13 30.81
        KGRBND(4) = 1                                                     40.13 30.81
        DO IXC = MXC, 1, -1                                               30.81
          IF (KGRPNT(IXC,1).GT.1) THEN                                    30.81
            KGRBND(3) = IXC                                               30.81
            GOTO 91                                                       30.81
          END IF                                                          30.81
  91    ENDDO                                                             30.81
        NGRBND = 2                                                        30.81
      ELSE                                                                30.81
        ALLOCATE(KGRPNTNEW(MXC,MYC))
        ALLOCATE(SCANND(MCGRD))
        NGRBND = 0
        DO ICGRD = 1, MCGRD
          SCANND(ICGRD) = 0
        ENDDO
        MCGRDOLD = MCGRD                                                  40.04
        DO 90 IXC = 1, MXC
          DO 80 IYC = 1, MYC
            ICGRD = KGRPNT(IXC,IYC)
            IF (ICGRD.GT.1) THEN
!             point with four wet neighbouring points which are all       40.04
!             no valid boundary points is eliminated.                     40.04
              IF ( NGRBGL.EQ.0 ) THEN                                     40.30
                IF ((.NOT. VALIDBP(IXC,IYC,KGRPNT,WNP)).OR.               40.04
     &             ((.NOT. VALIDBP(IXC-1,IYC,KGRPNT,WNP)).AND.            40.04
     &              (.NOT. VALIDBP(IXC+1,IYC,KGRPNT,WNP)).AND.            40.04
     &              (.NOT. VALIDBP(IXC,IYC-1,KGRPNT,WNP)).AND.            40.04
     &              (.NOT. VALIDBP(IXC,IYC+1,KGRPNT,WNP)))) THEN          40.04
                   KGRPNT(IXC,IYC) = 1                                    40.04
                   MCGRD = MCGRD - 1                                      40.04
                   SCANND(ICGRD) = 1                                      40.04
                   CALL MSGERR (1,                                        40.04
     &                          'Point removed from computational grid')  40.04
                   WRITE(PRINTF,17)  IXC, IYC                             40.04
  17               FORMAT ('Point with index (', I6, ',', I6,') was ',    40.04
     &                     'removed from computational grid, because ',   40.04
     &                     'it is part of a 1D configuration or 1D ',     40.04
     &                     'connection.')                                 40.04
                ENDIF                                                     40.04
              END IF                                                      40.30
              IF (SCANND(ICGRD).EQ.0) THEN
                IF (.NOT. PVALID(IXC-1,IYC,KGRPNT).AND.                   40.04
     &               VALIDBP(IXC,IYC,KGRPNT,WNP)) THEN                    40.04
                  IXOLD = IXC
                  IYOLD = IYC
                  NGRBND = NGRBND + 1
                  KGRBND(2*NGRBND-1) = IXOLD
                  KGRBND(2*NGRBND)   = IYOLD
                  SCANND(KGRPNT(IXOLD,IYOLD)) = 1
                  IDIR = 4                                                40.04
                  IF (WNP.EQ.3) THEN                                      40.04
                    CALL SEPARAREA(IXNEW, IYNEW, KGRPNT,IDIR)             40.04
                  ENDIF                                                   40.04
                  DO
                    IF (IDIR.EQ.1) THEN                                   40.04
                      IXNEW = IXOLD + 1                                   40.04
                      IYNEW = IYOLD                                       40.04
                    ELSE IF (IDIR.EQ.2) THEN                              40.04
                      IXNEW = IXOLD                                       40.04
                      IYNEW = IYOLD + 1                                   40.04
                    ELSE IF (IDIR.EQ.3) THEN                              40.04
                      IXNEW = IXOLD - 1                                   40.04
                      IYNEW = IYOLD                                       40.04
                    ELSE IF (IDIR.EQ.4) THEN                              40.04
                      IXNEW = IXOLD                                       40.04
                      IYNEW = IYOLD - 1                                   40.04
                    ENDIF                                                 40.04
                    IF (PVALID(IXNEW,IYNEW,KGRPNT) .AND.                  40.04
     &                  VALIDBP(IXNEW,IYNEW,KGRPNT,WNP) ) THEN            40.04
                      IF (SCANND(KGRPNT(IXNEW,IYNEW)) .GE. 1) GOTO 70
                      NGRBND = NGRBND + 1                                 40.04
                      KGRBND(2*NGRBND-1) = IXNEW                          40.04
                      KGRBND(2*NGRBND)   = IYNEW                          40.04
                      SCANND(KGRPNT(IXNEW,IYNEW)) = 1
                      IXOLD = IXNEW
                      IYOLD = IYNEW
                      IDIR  = IDIR - 1
                      IF (WNP.EQ.3) THEN                                  40.04
                        CALL SEPARAREA(IXNEW, IYNEW, KGRPNT,IDIR)         40.04
                      ENDIF                                               40.04
                      IF (IDIR.EQ.0) IDIR = 4                             40.04
                    ELSE
                      IDIR  = IDIR + 1
                      IF (.NOT. VALIDBP(IXNEW,IYNEW,KGRPNT,WNP).AND.      40.04
     &                    PVALID(IXNEW,IYNEW,KGRPNT) .AND.                40.30 40.04
     &                    NGRBGL.EQ.0 ) THEN                              40.30 40.04
                        KGRPNT(IXNEW,IYNEW) = 1                           40.04
                        MCGRD = MCGRD - 1                                 40.04
                        CALL MSGERR (1,                                   40.04
     &                          'Point removed from computational grid')  40.04
                        WRITE(PRINTF,18) IXNEW, IYNEW                     40.04
  18                    FORMAT ('Point with index (', I6 ,',',  I6,') ',  40.04
     &                          'was removed from computational grid ',   40.04
     &                          'because it is an isolated wet point ',   40.04
     &                          'or is part of a 1D configuration or ',   40.04
     &                          '1D connection.')                         40.04
                        SCANND(KGRPNT(IXNEW,IYNEW)) = 1                   40.04
                      END IF                                              40.04
                      IF (IDIR.EQ.5) IDIR = 1                             40.04
                    ENDIF
                  ENDDO
!                 the following indicates that curve is closed            40.04
  70              NGRBND = NGRBND + 1
                  KGRBND(2*NGRBND-1) = 0
                  KGRBND(2*NGRBND)   = 0
                ENDIF
              ENDIF
            ENDIF
  80      ENDDO
  90    ENDDO
        IF (MCGRD.NE.MCGRDOLD) THEN                                       40.04
          MCGRDNEW = 1                                                    40.04
          DO IXC = 1, MXC                                                 40.04
            DO IYC = 1, MYC                                               40.04
              IF (KGRPNT(IXC,IYC).NE.1) THEN                              40.04
                MCGRDNEW = MCGRDNEW + 1                                   40.04
                KGRPNTNEW(IXC,IYC) = MCGRDNEW                             40.04
              ELSE                                                        40.04
                KGRPNTNEW(IXC,IYC) = 1                                    40.04
              END IF                                                      40.04
            END DO                                                        40.04
          END DO                                                          40.04
          MCGRD = MCGRDNEW                                                40.04
          DO IXC = 1, MXC                                                 40.04
            DO IYC = 1, MYC                                               40.04
              KGRPNT(IXC,IYC) = KGRPNTNEW(IXC,IYC)                        40.04
            END DO                                                        40.04
          END DO                                                          40.04
        ENDIF                                                             40.04
        DEALLOCATE(KGRPNTNEW,SCANND)
      END IF                                                              30.81
      RETURN
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION PVALID (IX, IY, KGRPNT)
!                                                                      *
!***********************************************************************

      USE SWCOMM3                                                         40.41
!
      IMPLICIT NONE
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
!  1. Updates
!
!     New function for curvilinear version (ver. 40.00). May '98
!
!  2. Purpose
!
!     procedure to find whether a couple (ix,iy) represents
!     a valid grid point
!
!  3. Method
!
!     If one of the gridcounter (IX,IY) is less then 1 or greater than
!     maximum or point is an exception point the point is not valid.
!
!  4. Argument variables
!
!     IX, IY    input    x- and y-index of point under consideration
!     KGRPNT    input    indirect addresses for grid points
!
      INTEGER IX, IY, KGRPNT(MXC,MYC)
!
!  5. Parameter variables
!
!  6. Local variables
!
!     IENT    number of entries of this subroutine
!
      INTEGER ICGRD, IENT
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     CGBOUN
!     SEPARAREA
!     function VALIDBP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     PVALID = .TRUE.
!     If gridcounter in x-direction is less then 1 or greater than MXC or
!     If gridcounter in y-direction is less then 1 or greater than MYC or
!     If gridpoint is exception point PVALID = .FALSE.
!
! 13. Source text
!************************************************************************
!
      SAVE       IENT
      DATA       IENT /0/
      CALL STRACE (IENT, 'PVALID')
!
      PVALID = .TRUE.
      IF (IX.LT.1)   PVALID = .FALSE.
      IF (IY.LT.1)   PVALID = .FALSE.
      IF (IX.GT.MXC) PVALID = .FALSE.
      IF (IY.GT.MYC) PVALID = .FALSE.
      IF (PVALID) THEN
        ICGRD = KGRPNT(IX,IY)
        IF (ICGRD.LE.1) PVALID = .FALSE.
      ENDIF
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION VALIDBP (IX, IY, KGRPNT,WNP)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM3                                                         40.41
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
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
!     40.04  Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     August 2000 new function
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Check whether point with index (IX,IY) can be a valid boundary point
!
!  3. Method
!
!     The number of wet neighbouring points WNP is determined
!     Depending on this number certain configurations with wet and dry
!     points surrounding (IX,IY) are excluded
!
!  4. Argument variables
!
!     IX, IY    input    x- and y-index of point under consideration
!     KGRPNT    input    indirect addresses for grid points
!     WNP       output   number of wet neighbouring points
!
      INTEGER IX, IY, KGRPNT(MXC,MYC), WNP
!
!  5. Parameter variables
!
!  6. Local variables
!
!     IENT    number of entries of this subroutine
!
      INTEGER IENT
!
!  8. Subroutines used
!
!     Function PVALID
!
      LOGICAL PVALID
!
!  9. Subroutines calling
!
!     CGBOUN
!
! 10. Error messages
!
! 11. Remarks
!
!     This function prevends all projections of one cell width to
!     be a boundary point
!
! 12. Structure
!
!     Determine amount of Wet Neighbouring Points WNP
!     IF WNP =0  (isolated cell)     VALIDBP = .FALSE.
!     IF WNP =1  (one wet neighbour) VALIDBP = .FALSE.
!     IF WNP =2  (neighbouring points on straight line)
!                                    VALIDBP = .FALSE.
!
!      . W-d     (neighbouring points make angle with no wet point
!        | |      'between' them)    VALIDBP = .FALSE.
!      D-X-W
!        |
!      . D .     (X: point under consideration (assumed to be wet)
!                 W: wet neighbour
!                 D: dry neighbour  d: dry non-neighbour)
!                 .: either wet or dry point
!
!     IF WNP =3
!      . W-d     (point is 1D connection between areas or centre point
!        | |      of isolated 'half plus')
!      D-X-W                         VALIDBP = .FALSE.
!        | |
!      . W-d
!
!
! 13. Source text
!
!***********************************************************************
!
      SAVE       IENT
      DATA       IENT /0/
      CALL STRACE (IENT, 'VALIDBP')
!
      VALIDBP = .TRUE.
!      IF (PVALID(IX,IY,KGRPNT)) THEN
        WNP = 0
        IF (PVALID(IX-1,IY,KGRPNT)) WNP = WNP +1
        IF (PVALID(IX,IY-1,KGRPNT)) WNP = WNP +1
        IF (PVALID(IX+1,IY,KGRPNT)) WNP = WNP +1
        IF (PVALID(IX,IY+1,KGRPNT)) WNP = WNP +1
!         isolated point
        IF (WNP.EQ.0) VALIDBP = .FALSE.
!         point with one valid (wet) neighbouring grid point
        IF (WNP.EQ.1) VALIDBP = .FALSE.
!
!         neighbouring points on straight line (i.e. in fact 1D)
        IF ((WNP.EQ.2) .AND.(
     &      (PVALID(IX,IY-1,KGRPNT).AND.PVALID(IX,IY+1,KGRPNT)) .OR.
     &      (PVALID(IX-1,IY,KGRPNT).AND.PVALID(IX+1,IY,KGRPNT)) .OR.
!         neighbouring points make angle but no wet point 'between' them
     &      (PVALID(IX-1,IY,KGRPNT).AND.PVALID(IX,IY+1,KGRPNT) .AND.
     &       .NOT. PVALID(IX-1,IY+1,KGRPNT)) .OR.
     &      (PVALID(IX-1,IY,KGRPNT).AND.PVALID(IX,IY-1,KGRPNT) .AND.
     &       .NOT. PVALID(IX-1,IY-1,KGRPNT)) .OR.
     &      (PVALID(IX+1,IY,KGRPNT).AND.PVALID(IX,IY-1,KGRPNT) .AND.
     &       .NOT. PVALID(IX+1,IY-1,KGRPNT)) .OR.
     &      (PVALID(IX+1,IY,KGRPNT).AND.PVALID(IX,IY+1,KGRPNT) .AND.
     &       .NOT. PVALID(IX+1,IY+1,KGRPNT)) )  )  VALIDBP = .FALSE.
!
!        point (IX,IY) is 1D connection between areas or isolated
!        centre point of 'half plus'
        IF ((WNP.EQ.3) .AND.(
     &       (.NOT.PVALID(IX-1,IY,KGRPNT) .AND.
     &        .NOT.PVALID(IX+1,IY-1,KGRPNT).AND.
     &        .NOT.PVALID(IX+1,IY+1,KGRPNT)).OR.
     &       (.NOT.PVALID(IX+1,IY,KGRPNT) .AND.
     &        .NOT.PVALID(IX-1,IY-1,KGRPNT).AND.
     &        .NOT.PVALID(IX-1,IY+1,KGRPNT)).OR.
     &       (.NOT.PVALID(IX,IY-1,KGRPNT) .AND.
     &        .NOT.PVALID(IX-1,IY+1,KGRPNT).AND.
     &        .NOT.PVALID(IX+1,IY+1,KGRPNT)).OR.
     &       (.NOT.PVALID(IX,IY+1,KGRPNT) .AND.
     &        .NOT.PVALID(IX-1,IY-1,KGRPNT).AND.
     &        .NOT.PVALID(IX+1,IY-1,KGRPNT)) )  )  VALIDBP = .FALSE.
!
      RETURN
      END
!
!*******************************************************************
!                                                                  *
      SUBROUTINE SEPARAREA(IX, IY, KGRPNT,IDIR)
!                                                                  *
!*******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
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
!     40.04  Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     August 2000 new subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Separate areas that could be connected with a one cell connection
!
!  3. Method
!
!     Redirect original IDIR ( '<#' in plot below) to new direction('=>')
!
!     .    D    W -- W         X: point under consideration (wet)
!               |    |         W: wet point
!     W -- W <# X => W         D: dry point
!     |    |    |              .: either wet or dry point
!     W -- W    D    .
!
!     for all 8 different situations (4 rotation variations and their
!                                     mirrored situations)
!
!  4. Argument variables
!
!     IX, IY    input          x- and y-index of point under consideration
!     KGRPNT    input          indirect addresses for grid points
!     IDIR      input/output   index for direction (see subroutine CGBOUN)
!
      INTEGER IX, IY, KGRPNT(MXC,MYC), IDIR
!
!  5. Parameter variables
!
!  6. Local variables
!
!     IENT    number of entries of this subroutine
!
      INTEGER IENT
!
!  8. Subroutines used
!
!     Function PVALID
!
      LOGICAL PVALID
!
!  9. Subroutines calling
!
!     CGBOUN
!
! 10. Error messages
!
! 11. Remarks
!
!     This subroutine is only used if number of Wet Neighbouring
!     Points WNP = 3
!
! 12. Structure
!
!     In this plot the first situation is illustrated
!
!     .    d    W -- W         X  : point under consideration (wet)
!               |    |         W  : wet point
!     W -- W <# X => W         D,d: dry point
!     |    |    |              .  : either wet or dry point
!     W -- W    D    .
!
!
!     IF dry neighbour point D is below X and upper left point is dry
!     redirect IDIR to 1
!     IF dry neighbour point D is right of X and lower left point is dry
!     redirect IDIR to 2
!     IF dry neighbour point D is above X and lower right point is dry
!     redirect IDIR to 3
!     IF dry neighbour point D is left of X and upper right point is dry
!     redirect IDIR to 4
!
!     mirrored situation
!
!     IF dry neighbour point D is above X and lower left point is dry
!     redirect IDIR to 4
!     IF dry neighbour point D is left of X and lower right point is dry
!     redirect IDIR to 1
!     IF dry neighbour point D is below X and upper right point is dry
!     redirect IDIR to 2
!     IF dry neighbour point D is right of X and upper left point is dry
!     redirect IDIR to 3
!
! 13. Source text
!
!***********************************************************************
!
      DATA      IENT /0/
      CALL STRACE (IENT,'SEPARAREA')
!
!     In case there are 3 wet neighbouring points and validbp(ix,iy,.)
!
      IF(.NOT.PVALID(IX-1,IY+1,KGRPNT) .AND.
     &    .NOT.PVALID(IX,IY-1,KGRPNT)) IDIR = 1
      IF(.NOT.PVALID(IX-1,IY-1,KGRPNT) .AND.
     &    .NOT.PVALID(IX+1,IY,KGRPNT)) IDIR = 2
      IF(.NOT.PVALID(IX+1,IY-1,KGRPNT) .AND.
     &    .NOT.PVALID(IX,IY+1,KGRPNT)) IDIR = 3
      IF(.NOT.PVALID(IX+1,IY+1,KGRPNT) .AND.
     &    .NOT.PVALID(IX-1,IY,KGRPNT)) IDIR = 4
!
!     mirrored situation
      IF(.NOT.PVALID(IX-1,IY-1,KGRPNT) .AND.
     &    .NOT.PVALID(IX,IY+1,KGRPNT)) IDIR = 4
      IF(.NOT.PVALID(IX+1,IY-1,KGRPNT) .AND.
     &    .NOT.PVALID(IX-1,IY,KGRPNT)) IDIR = 1
      IF(.NOT.PVALID(IX+1,IY+1,KGRPNT) .AND.
     &    .NOT.PVALID(IX,IY-1,KGRPNT)) IDIR = 2
      IF(.NOT.PVALID(IX-1,IY+1,KGRPNT) .AND.
     &    .NOT.PVALID(IX+1,IY,KGRPNT)) IDIR = 3
!
      RETURN
      END
!
!*******************************************************************
!                                                                  *
      SUBROUTINE INITVA( AC2, SPCSIG, SPCDIR, KGRPNT )                    40.80 40.31
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE TIMECOMM                                                        40.41
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
!PUN      USE SIZES                                                           40.95
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
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.03, 40.13: Nico Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.62: Ben Payment (MSU), Tim Campbell (NRL)
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.70, Oct. 97: New subroutine
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of SPCDIR
!     30.82, Dec. 98: Corrected the arguments in CALL SINTRP(..)
!     34.01, Feb. 99: Introducing STPNOW
!     40.00, Aug. 99: modification for 1D mode; new option INIT PAR
!                     init restart added
!     40.03, Nov. 99: after reading comment line, jump to 110 (not 100)
!                     additional test output added
!                     possibility added to initialize in limited region (PAR case)
!                     function EQCSTR used to compare strings
!     40.13, Jan. 01: option Spherical was not yet taken care of
!                     ! is now allowed as comment sign in a restart file
!     40.13, Oct. 01: error message removed, command MODE not required any more
!     40.31, Oct. 03: small changes
!     40.31, Dec. 03: appending number to file name i.c. of
!                     parallel computing
!     40.41, Sep. 04: small changes
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.62, Jul. 06: modified HOTSTART functionality to handle reading from a
!                     single hotfile or from multiple hotfiles when running in
!                     parallel with MPI
!     40.80, Jun. 07: extension to unstructured grids
!
!  2. Purpose
!
!     process command INIT and compute initial state of the wave field
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!  8. Subroutines used
!
      LOGICAL :: STPNOW, EQCSTR                                           40.03
!
! 13. Source text
!
!
      INTEGER    KGRPNT(MXC,MYC)                                          40.00
      INTEGER JXMAX, JYMAX, NPTOT                                         40.80 40.62
      REAL       AC2(MDC,MSC,MCGRD)                                       40.31
      LOGICAL SINGLEHOT, PTNSUBGRD                                        40.62
      LOGICAL    KEYWIS, LERR                                             40.31 40.00
      CHARACTER  RLINE *80                                                40.00
      INTEGER IID, IUNITAC                                                40.41
      REAL    ACTMP(MDC), DIRTMP(MDC)                                     40.41
      LOGICAL EQREAL                                                      40.41
!
      SAVE       IENT
      DATA       IENT /0/
      CALL STRACE (IENT, 'INITVA')                                        40.00
!
!     ------------------------------------------------------------------
!
!     ***'initial conditions'  Definition of initial conditions  ***
!     *** for MODE DYNAMIC
!
! ============================================================
!
!               | -> DEFault
!               |
!     INITial  <   ZERO                                                   40.00
!               |
!               |  PAR  [hs] [per] [dir] [dd]
!               |
!               |            | -> MULTiple |                              40.62
!               |  HOTStart <               >  'fname'                    40.62 40.00
!               |            |    SINGle   |                              40.62
!
! ============================================================
!
      CALL IGNORE ('COND')
      LERR =  .FALSE.                                                     40.00
!     error message removed because MODE is not required any more,        40.13
!     and INIT can be useful also in stationary mode                      40.13
      CALL INKEYW ('STA', 'DEF')
      IF (KEYWIS('PAR')) THEN
!       initial state defined by wave parameters
        IF (MXC.LE.0 .AND. OPTG.NE.5) THEN                                40.80
          CALL MSGERR (2,
     &        'command INIT should follow CGRID')                         40.80
          LERR = .TRUE.
        ENDIF
        IF (MCGRD .LE. 1 .AND. nverts.LE.0) THEN                          40.80
          CALL MSGERR (2,
     &        'command INIT should follow READ BOT or READ UNSTRUC')      40.80
          LERR = .TRUE.
        ENDIF
        ICOND = 2
        CALL INREAL ('HSIG', SPPARM(1), 'REQ', 0.)
        CALL INKEYW ('STA', '    ')                                       40.00
        IF (KEYWIS('MEAN')) THEN
          IF (FSHAPE.GT.0) FSHAPE=-FSHAPE
        ELSE IF (KEYWIS('PEAK')) THEN                                     40.00
          IF (FSHAPE.LT.0) FSHAPE=-FSHAPE                                 40.00
        ENDIF
        CALL INREAL('PER', SPPARM(2), 'REQ', 0.)
        IF (SPPARM(2).LE.0.) CALL MSGERR (2, 'Period must be >0')
        IF (SPPARM(2).GT.PI2/SPCSIG(1)) THEN
          CALL MSGERR (2,
     &          'Inc. freq. lower than lowest spectral freq.')
          WRITE(PRINTF,11) 1./SPPARM(2),' < ',SPCSIG(1)/PI2,
     &                     ' lowest'
        ENDIF
        IF (SPPARM(2).LT.PI2/SPCSIG(MSC)) THEN
          CALL MSGERR (2,
     &          'Inc. freq. higher than highest spectral freq.')
          WRITE(PRINTF,11) 1./SPPARM(2),' > ',SPCSIG(MSC)/PI2,
     &                     'highest'
        ENDIF
 11     FORMAT(' Inc. freq. = ',F9.5, A3, F9.5, '=',A7,' freq')
        CALL INREAL('DIR',  SPPARM(3), 'REQ', 0.)
        IF (DSHAPE.EQ.1) THEN
          CALL INREAL('DD', SPPARM(4), 'STA', 30.)
        ELSE
          CALL INREAL('DD', SPPARM(4), 'STA', 2.)
        ENDIF
!       give boundaries of region where the initial condition applies     40.03
        IF (OPTG.NE.5) THEN                                               40.80
           CALL ININTG ('IX1', IX1, 'STA', 0)                             40.03
           CALL ININTG ('IX2', IX2, 'STA', MXC-1)                         40.03
           CALL ININTG ('IY1', IY1, 'STA', 0)                             40.03
           CALL ININTG ('IY2', IY2, 'STA', MYC-1)                         40.03
        ENDIF                                                             40.80
        IF (.NOT.LERR) THEN
          CALL SSHAPE (AC2(1,1,1), SPCSIG, SPCDIR, FSHAPE, DSHAPE)
!         copy computed spectrum to all internal grid points
          IF (OPTG.NE.5) THEN                                             40.80
!            --- structured grid                                          40.80
             DO IX = IX1+1, IX2+1                                         40.03
               DO IY = IY1+1, IY2+1                                       40.03
                 INDX = KGRPNT(IX,IY)
                 IF (INDX.GT.1)
     &           CALL SINTRP (1., 0., AC2(1,1,1), AC2(1,1,1),
     &                          AC2(1,1,INDX),SPCDIR,SPCSIG)              30.82
               ENDDO                                                      40.03
             ENDDO                                                        40.00
!            reset action density AC2(*,*,1) to 0
             DO ID = 1, MDC
               DO IS = 1, MSC
                 AC2(ID,IS,1) = 0.
               ENDDO
             ENDDO
          ELSE                                                            40.80
!            --- unstructured grid                                        40.80
             DO K = 2, nverts                                             40.80
                CALL SINTRP (1., 0., AC2(1,1,1), AC2(1,1,1),              40.80
     &                       AC2(1,1,K),SPCDIR,SPCSIG)                    40.80
             ENDDO                                                        40.80
          ENDIF                                                           40.80
        ENDIF
      ELSE IF (KEYWIS('ZERO')) THEN
!
!       zero initial state
!
!       the statements below work also for unstructured grids             40.80
!       since, MCGRD = nverts (see SwanInitCompGrid)                      40.80
        DO INDX = 1, MCGRD                                                40.00
          DO ID = 1, MDC
            DO IS = 1, MSC
              AC2(ID,IS,INDX) = 0.
            ENDDO
          ENDDO
        ENDDO
        ICOND = 3
      ELSE IF (KEYWIS('HOTS') .OR. KEYWIS('REST')) THEN                   40.00
        IID       = 0                                                     40.41
        IUNITAC   = 0                                                     40.41
        ACTMP     = 0.                                                    40.41
        DIRTMP(:) = SPCDIR(:,1)                                           40.41
!       initialize using spectra from a HOTFILE                           40.00
        IF (MXC.LE.0 .AND. OPTG.NE.5) CALL MSGERR (2,                     40.80
     &        'command INIT should follow CGRID')                         40.80
        IF (MCGRD.LE.1 .AND. nverts.LE.0) CALL MSGERR (2,                 40.80
     &        'command INIT should follow READ BOT or READ UNSTRUC')      40.80
        ICOND = 4
        CALL INKEYW ('STA', 'MULT')                                       40.62
        IF (KEYWIS ('SING')) THEN                                         40.62
          SINGLEHOT = .TRUE.                                              40.62
          JXMAX = MXCGL                                                   40.62
          JYMAX = MYCGL                                                   40.62
          CALL IGNORE ('SING')                                            40.62
        ELSE                                                              40.62
          SINGLEHOT = .FALSE.                                             40.62
          JXMAX = MXC                                                     40.62
          JYMAX = MYC                                                     40.62
          CALL IGNORE ('MULT')                                            40.62
        END IF                                                            40.62
        NPTOT = JXMAX*JYMAX                                               40.80
        IF (OPTG.EQ.5) NPTOT = nverts                                     40.80
        CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
!       --- append node number to FILENM in case of parallel computing    40.31
        IF ( PARLL .AND. .NOT.SINGLEHOT ) THEN                            40.62 40.31
           ILPOS = INDEX ( FILENM, ' ' )-1                                40.31
           WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                        40.31
  33       FORMAT('-',I3.3)                                               40.31
        END IF                                                            40.31
!PUN        IF (.NOT.SINGLEHOT) FILENM= TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)  40.95
        NREF   = 0
        IOSTAT = 0
        CALL FOR (NREF, FILENM, 'OF', IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
 100    READ (NREF, 102) RLINE
 102    FORMAT (A)                                                        40.00
        IF (RLINE(1:4).NE.'SWAN') CALL MSGERR (3, FILENM//
     &        ' is not a correct hotstart file')
 110    READ (NREF, 102) RLINE
        IF (RLINE(1:1).EQ.COMID .OR. RLINE(1:1).EQ.'!') GOTO 110          40.13
        IF (EQCSTR(RLINE,'TIME')) THEN
          READ (NREF, *) IIOPT
          IF (ITEST.GE.50) WRITE (PRTEST, 122) IIOPT                      40.03
 122      FORMAT (' time coding option:', I2)
          READ (NREF, 102) RLINE
!         in stationary mode, warning
          IF (NSTATM.EQ.0) CALL MSGERR (1,
     &                  'Time info in hotfile ignored')                   40.03
        ELSE
          IIOPT = -1
          IF (NSTATM.EQ.1) CALL MSGERR (1,
     &                  'No time info in hotfile')                        40.03
        ENDIF
        IF (EQCSTR(RLINE,'LOCA') .OR. EQCSTR(RLINE,'LONLAT')) THEN        40.13
          READ (NREF, *) NUMPTS
          IF (NUMPTS.NE.NPTOT) THEN                                       40.80 40.62
            CALL MSGERR (2,
     &      'grid on hotstart file differs from one in CGRID command')    40.00
            WRITE (PRINTF, 123) NPTOT, NUMPTS                             40.80 40.62
 123        FORMAT (1X, I8, ' points in comp.grid; on file:', I8)         40.03
          ENDIF
          IF (ITEST.GE.50) WRITE (PRTEST, 124) NUMPTS                     40.03
 124      FORMAT (1X, I8, '  output locations')
          DO IP = 1, NUMPTS
            READ (NREF, *)
          ENDDO
          READ (NREF, 102) RLINE
        ENDIF
        IF (EQCSTR(RLINE(2:5),'FREQ')) THEN                               40.03
          READ (NREF, *) NUMFRE
          IF (NUMFRE.NE.MSC) CALL MSGERR (2,
     &    'grid on hotstart file differs from one in CGRID command')      40.00
          IF (ITEST.GE.50) WRITE (PRTEST, 126) NUMFRE                     40.03
 126      FORMAT (1X, I6, '  frequencies')
          DO IP = 1, NUMFRE
            READ (NREF, *)
          ENDDO
          READ (NREF, 102) RLINE
        ENDIF
        IF (EQCSTR(RLINE(2:4),'DIR')) THEN                                40.03
          READ (NREF, *) NUMDIR
          IF (NUMDIR.NE.MDC) CALL MSGERR (2,
     &    'grid on hotstart file differs from one in CGRID command')      40.00
          IF (ITEST.GE.50) WRITE (PRTEST, 128) NUMDIR                     40.03
 128      FORMAT (1X, I6, '  directions')
          DO IP = 1, NUMDIR
            READ (NREF, *) DIRTMP(IP)                                     40.41
          ENDDO
          IF (.NOT.EQREAL(DIRTMP(1)*PI/180.,SPCDIR(1,1)))                 40.41
     &       IID = NINT(REAL(MDC)*(DIRTMP(1)-ALPC)/360.)                  40.41
          READ (NREF, 102) RLINE
        ENDIF
        READ (NREF, *) NQUA
        IF (NQUA.NE.1) CALL MSGERR (2,'NQUA>1: incorrect hotstart file')  40.00
        READ (NREF, 102) RLINE
        IF (ITEST.GE.50) WRITE (PRTEST, 130) RLINE                        40.03
 130    FORMAT (1X, 'quantity: ', A)
        READ (NREF, 102) RLINE                                            40.00
        IF (EQCSTR(RLINE(3:3),'S')) IUNITAC = 1                           40.41
        READ (NREF, 102) RLINE                                            40.00
!
!       reading of heading is completed, read time if nonstationary
!
        IF (IIOPT.GE.0) THEN
          READ (NREF, 102) RLINE                                          40.00
          CALL DTRETI (RLINE(1:18), IIOPT, TIMCO)                         40.00
          WRITE (PRINTF, 210) RLINE(1:18)
 210      FORMAT (' initial condition read for time: ', A)
        ENDIF
!
        IF (OPTG.NE.5) THEN                                               40.80
!
!       --- structured grid                                               40.80
!
          DO 290 JX = 1, JXMAX                                            40.62
            IF (SINGLEHOT) THEN                                           40.62
               IX = JX-MXF+1
            ELSE
               IX = JX
            ENDIF
            DO 280 JY = 1, JYMAX                                          40.62
              IF (SINGLEHOT) THEN                                         40.62
                 IY = JY-MYF+1
              ELSE
                 IY = JY
              ENDIF
              PTNSUBGRD = .TRUE.                                          40.62
              IF ( SINGLEHOT .AND.                                        40.62
     &          (MXF.GT.JX .OR. MXL.LT.JX .OR.MYF.GT.JY .OR. MYL.LT.JY))  40.62
     &           PTNSUBGRD = .FALSE.                                      40.62
!
              READ (NREF, 102) RLINE
              INDX = KGRPNT(IX,IY)
              IF (INDX.EQ.1 .AND. PTNSUBGRD) THEN                         40.62
                IF (RLINE(1:6).NE.'NODATA') THEN
                   CALL MSGERR (2,
     &                'valid spectrum for non-existing grid point')
                   WRITE (PRINTF, *) IX-1, IY-1
                ENDIF
              ELSE
                IF (EQCSTR(RLINE,'NODATA').OR.EQCSTR(RLINE,'ZERO')) THEN
                  IF (PTNSUBGRD) THEN                                     40.62
                     DO IS = 1, MSC
                       DO ID = 1, MDC
                          AC2(ID,IS,INDX) = 0.
                       ENDDO
                     ENDDO
                     IF (ITEST.GE.150) WRITE (PRTEST, 222) IX-1, IY-1     40.03
                  ENDIF
 222              FORMAT (' zero spectrum or no data for point:', 2I4)    40.03
                ELSE
!                 first determine factor
                  READ (NREF, *) AFAC
!                 multiply with factor to account for transition from     40.00
!                 energy/Hz/degr to energy/(2*pi rad/s)/rad
                  AFAC = AFAC * 90. / (PI**2)
                  DO IS = 1, MSC
                    AFAC1 = AFAC/SPCSIG(IS)                               40.41
!                   Read spectral energy densities from file
                    READ (NREF, *) (ACTMP(ID), ID=1,MDC)                  40.41
                    IF (.NOT.PTNSUBGRD) CYCLE                             40.62
                    IF (IUNITAC.EQ.1) THEN                                40.41
                       DO ID = 1, MDC
                         J = MODULO ( IID - 1 + ID , MDC ) + 1            40.69 40.41
                         AC2(J,IS,INDX) = AFAC * ACTMP(ID)                40.41
                       ENDDO
                    ELSE
                       DO ID = 1, MDC
                         J = MODULO ( IID - 1 + ID , MDC ) + 1            40.69 40.41
                         AC2(J,IS,INDX) = AFAC1 * ACTMP(ID)               40.41
                       ENDDO
                    END IF
                  ENDDO
                  IF (ITEST.GE.150) WRITE (PRTEST, 224) IX-1, IY-1, AFAC  40.03
 224              FORMAT (' spectrum in point:', 2I4,'  factor=', E12.4)  40.03
                ENDIF
              ENDIF
 280        CONTINUE
 290      CONTINUE
        ELSE                                                              40.80
!
!          --- unstructured grids                                         40.80
!
           DO K = 1, nverts
              READ (NREF, 102) RLINE
              IF (EQCSTR(RLINE,'NODATA') .OR. EQCSTR(RLINE,'ZERO')) THEN
                 DO IS = 1, MSC
                    DO ID = 1, MDC
                       AC2(ID,IS,K) = 0.
                    ENDDO
                 ENDDO
                 IF (ITEST.GE.150) WRITE (PRTEST, 223) K
 223             FORMAT (' zero spectrum or no data for vertex:', I6)
              ELSE
!                first determine factor
                 READ (NREF, *) AFAC
!                multiply with factor to account for transition from
!                energy/Hz/degr to energy/(rad/s)/rad
                 AFAC = AFAC * 90. / (PI**2)
                 DO IS = 1, MSC
                    AFAC1 = AFAC/SPCSIG(IS)
!                   Read spectral energy densities from file
                    READ (NREF, *) (ACTMP(ID), ID=1,MDC)
                    IF (IUNITAC.EQ.1) THEN
                       DO ID = 1, MDC
                         J = MOD ( IID - 1 + ID , MDC ) + 1
                         AC2(J,IS,K) = AFAC * ACTMP(ID)
                       ENDDO
                    ELSE
                       DO ID = 1, MDC
                         J = MOD ( IID - 1 + ID , MDC ) + 1
                         AC2(J,IS,K) = AFAC1 * ACTMP(ID)
                       ENDDO
                    END IF
                 ENDDO
                 IF (ITEST.GE.150) WRITE (PRTEST, 225) K, AFAC
 225             FORMAT (' spectrum in vertex:', I6, '  factor=', E12.4)
              ENDIF
           ENDDO
        ENDIF                                                             40.80
        CLOSE (NREF)
      ELSE
!       default initial wave state, will be computed later by subr SWINCO
        CALL IGNORE ('DEF')
        ICOND = 1
      ENDIF
      RETURN
!     end of subr INITVA
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE BACKUP (AC2, SPCSIG, SPCDIR, KGRPNT,
     &                   XCGRID, YCGRID)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
!PUN      USE SIZES                                                           40.95
!ADC      USE Couple2Swan, ONLY: SwanHotStartUnit                             41.20
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
!     ver 40.00, Apr 1998 by N.Booij: new subroutine
!
!   Purpose
!
!  0. Authors
!
!     30.82: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     34.01: Jeroen Adema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!     41.20: Casey Dietrich
!
!  1. Updates
!
!     40.00, Apr. 98: New subroutine
!     30.82, Oct. 98: Updated description of SPCDIR
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Nov. 99: ITMOPT is written as time coding option
!     40.13, Jan. 01: option Spherical was not yet taken care of
!     40.31, Dec. 03: appending number to file name i.c. of
!                     parallel computing
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Jun. 07: extension to unstructured grids
!     41.20, Mar. 10: extension to tightly coupled ADCIRC+SWAN model
!
!  2. Purpose
!
!     backup current state of the wave field to a file
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
! i   XCGRID: Coordinates of computational grid in x-direction            30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.82
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.82
!
!  8. Subroutines used
!
      LOGICAL STPNOW                                                      34.01
!
! 13. Source text
!
      REAL     AC2(MDC,MSC,MCGRD)
      INTEGER  KGRPNT(MXC,MYC)
      CHARACTER (LEN=8) :: CRFORM = '(2F14.4)'                            40.41
      LOGICAL  EQREAL                                                     40.41
      SAVE     IENT
      DATA     IENT /0/
      CALL STRACE (IENT, 'BACKUP')
!
!     ==================================================================
!
!     HOTFile  'FNAME'                                                    40.00
!
!     ==================================================================
!
      CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
!ADC      IF ( SwanHotStartUnit.EQ.67 ) THEN                                  41.20
!ADC         FILENM = "swan.67"
!ADC         SwanHotStartUnit = 68
!ADC      ELSE
!ADC         FILENM = "swan.68"
!ADC         SwanHotStartUnit = 67
!ADC      ENDIF
!     --- append node number to FILENM in case of parallel computing      40.31
      IF ( PARLL ) THEN                                                   40.31
         ILPOS = INDEX ( FILENM, ' ' )-1                                  40.31
         WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                          40.31
  33     FORMAT('-',I3.3)                                                 40.31
      END IF                                                              40.31
!PUN      FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                       40.95
      NREF   = 0
      IOSTAT = 0
      CALL FOR (NREF, FILENM, 'UF', IOSTAT)
      IF (STPNOW()) RETURN                                                34.01
      WRITE (NREF, 102) 'SWAN   1', 'Swan standard file, version'
      IF (NSTATM.EQ.1) THEN
        WRITE (NREF, 102) 'TIME', 'time-dependent data'
 102    FORMAT (A, T41, A)                                                40.00
        WRITE (NREF, 103) ITMOPT, 'time coding option'                    40.03
 103    FORMAT (I6, T41, A)                                               40.00
      ENDIF
      IF (KSPHER.EQ.0) THEN                                               40.13
        WRITE (NREF, 102) 'LOCATIONS', 'locations in x-y-space'
        CRFORM = '(2F14.4)'                                               40.41
      ELSE                                                                40.13
        WRITE (NREF, 102) 'LONLAT', 'locations on the globe'              40.13
        CRFORM = '(2F12.6)'                                               40.41
      ENDIF                                                               40.13
      IF (OPTG.NE.5) THEN                                                 40.80
         WRITE (NREF, 104) MXC*MYC, MXC, MYC, 'number of locations'
 104     FORMAT (I8, 2I6, T41, A)
         DO IX = 1, MXC
           DO IY = 1, MYC
             IF ( EQREAL(XCGRID(IX,IY), EXCFLD(8)) .AND.                  40.41
     &            EQREAL(YCGRID(IX,IY), EXCFLD(9)) ) THEN                 40.41
               WRITE (NREF, FMT=CRFORM) DBLE(EXCFLD(8)) + DBLE(XOFFS),
     &                                  DBLE(EXCFLD(9)) + DBLE(YOFFS)
             ELSE
             WRITE (NREF, FMT=CRFORM) DBLE(XCGRID(IX,IY)) + DBLE(XOFFS),
     &                                DBLE(YCGRID(IX,IY)) + DBLE(YOFFS)   40.00
             ENDIF                                                        40.41
           ENDDO
         ENDDO
      ELSE                                                                40.80
         WRITE (NREF, 105) nverts, 'number of locations'                  40.80
 105     FORMAT (I8, T41, A)
         DO K = 1, nverts                                                 40.80
            WRITE (NREF, FMT=CRFORM) DBLE(xcugrd(K)) + DBLE(XOFFS),       40.80
     &                               DBLE(ycugrd(K)) + DBLE(YOFFS)        40.80
         ENDDO                                                            40.80
      ENDIF                                                               40.80
      WRITE (NREF, 102) 'RFREQ', 'relative frequencies in Hz'             40.00
      WRITE (NREF, 103) MSC, 'number of frequencies'                      40.00
      DO 120 IS = 1, MSC
        WRITE (NREF, 114) SPCSIG(IS)/PI2
 114    FORMAT (F10.4)
 120  CONTINUE
      WRITE (NREF, 102) 'CDIR', 'spectral Cartesian directions in degr'   40.00
      WRITE (NREF, 103) MDC, 'number of directions'
      DO 130 ID = 1, MDC
        WRITE (NREF, 124) SPCDIR(ID,1)*180./PI                            30.82
 124    FORMAT (F10.4)
 130  CONTINUE
      WRITE (NREF, 132) 1
 132  FORMAT ('QUANT', /, I6, T41, 'number of quantities in table')       40.00
      WRITE (NREF, 102) 'AcDens', 'action densities'
      WRITE (NREF, 102) 'm2s/Hz/deg', 'unit'                              40.31 40.00
      WRITE (NREF, 102) '0.',     'exception value'                       40.00
!
!     writing of heading is completed, write time if nonstationary
!
      IF (NSTATM.EQ.1) THEN
        WRITE (NREF, 202) CHTIME                                          40.00
 202    FORMAT (A18, T41, 'date and time')
      ENDIF
!
      IF (OPTG.NE.5) THEN                                                 40.80
         DO 290 IX = 1, MXC
           DO 280 IY = 1, MYC
             INDX = KGRPNT(IX,IY)
             IF (INDX.EQ.1) THEN
               WRITE (NREF,220) 'NODATA'                                  40.08
             ELSE
               CALL WRSPEC (NREF, AC2(1,1,INDX))                          40.00
             ENDIF
 280       CONTINUE
 290     CONTINUE
      ELSE                                                                40.80
         DO K = 1, nverts                                                 40.80
            CALL WRSPEC (NREF, AC2(1,1,K))                                40.80
         ENDDO                                                            40.80
      ENDIF                                                               40.80
 220  FORMAT (A6)                                                         40.08
      CLOSE (NREF)
      RETURN
!     end of subr BACKUP
      END
