!
!     SWAN/SWREAD  file 2 of 2
!
!  Contents of this file:
!     SPROUT: Reading and processing of the user output commands
!     SWREPS: Reading and processing of the commands defining output points
!     SWREOQ: Reading and processing of the output requests
!     SIRAY : Searching the first point on a ray where the depth is DP
!     SWNMPS
!     SVARTP
!     SWBOUN                                                              40.00
!     BCFILE                                                              40.00
!     BCWAMN                                                              40.00
!     BCWW3N
!     SWBCPT
!     RETSTP                                                              40.00
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SPROUT (FOUND, BOTLEV, WATLEV)                           40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
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
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     32.02: Roeland Ris & Cor van der Schelde (1D version)
!     34.01: Jeroen Adema
!     40.02: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!    100.04, Nov. 92: Filename of plotfile will be given by user
!     30.70, Nov. 97: Arguments BOTLEV and WATLEV added
!     32.02, Feb. 98: 1D version introduced
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Apr. 98: Removed reference to commons KAART and KAR
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.81, Nov. 98: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     30.81, Jan. 99: Replaced variable FROM by FROM_ (because FROM is a
!                     reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Sep. 00: inconsistency with manual corrected
!     40.02, Oct. 00: Initialisation of IERR
!     40.31, Nov. 03: removing POOL construction and HPGL functionality
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Reading and processing of the user output commands
!
!  3. Method
!
!     If the first characters of the last read command are equal to a
!     given string (KEYWIS ('STRING')), the keywords and varia-
!     bles of this command are further read and processed
!
!  4. Argument variables
!
! i   BOTLEV: Bottom levels                                               30.70
! i   WATLEV: Water levels                                                30.70
!
      REAL    BOTLEV(*)                                                   30.70
      REAL    WATLEV(*)                                                   30.70
!
!  6. Local variables
!
!  8. Subroutines used
!
!     MSGERR
!     SWREPS
!     SWREOQ
!     STPNOW
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
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!             Most of the source code will be clear with the
!             aid of the user manual, the system documentation
!             and the additional comments in the source code.
!     ----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL   FOUND
      LOGICAL   KEYWIS
      SAVE IENT
      DATA IENT/0/                                                        40.31 30.81
      CALL STRACE (IENT,'SPROUT')
!
      FOUND = .FALSE.
!
!     definition of output point sets
!
      CALL SWREPS ( FOUND, BOTLEV, WATLEV )                               40.31
      IF (STPNOW()) RETURN                                                34.01
      IF (FOUND) RETURN
!
!     output requests
!
      CALL SWREOQ ( FOUND )                                               40.31 30.90
      IF (STPNOW()) RETURN                                                34.01
      IF (FOUND) RETURN
!
      IF (KEYWIS('SIT') .OR. KEYWIS('PLA')) THEN
        CALL MSGERR(2,'Keyword SITES is no longer maintained')            40.31
        GOTO 800                                                          40.31
      ENDIF
!
      IF (KEYWIS ('LIN')) THEN
        CALL MSGERR(2,'Keyword LINE is no longer maintained')             40.31
        GOTO 800                                                          40.31
      ENDIF
!     -------------------------------------------------------------------
!     ***** command name not found *****
      RETURN
!
 800  FOUND = .TRUE.
      RETURN
! *   end of subroutine SPROUT *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWREPS ( FOUND, BOTLEV, WATLEV )                         40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.31
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
!     30.70, 40.03, 40.13: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde (1D version)
!     34.01: Jeroen Adema
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.70, Nov. 97: comm ISO, inquire pointer added to get correct value
!                     for IADRAY
!     30.70, Nov. 97: comm ISO, offset origin added in message concerning rays
!                     declaration INT SIRAY added
!     30.70, Nov. 97: arguments BOTLEV and WATLEV added
!     30.72, Feb. 98: Declaration of Argument variables updated
!     32.02, Feb. 98: 1D version introduced
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Apr. 98: removed reference to commons KAART and KAR
!     30.81, Nov. 98: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.01, Sep. 99: XASM and YASM replace fixed numbers
!     33.09, Sep. 00: modifications in view of spherical coordinates
!     40.03, Sep. 00: inconsistency with manual corrected
!     40.13, Sep. 01: nesting in curvilinear grid: division by 0 prevented
!     40.30, May  03: introduction distributed-memory approach using MPI
!     40.31, Dec. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Mar. 08: extension to unstructured grids
!
!  2. PURPOSE
!
!     Reading and processing of the commands defining output points
!
!  4. Argument variables (updated 30.72)
!
!     BOTLEV: input  bottom levels                                        30.70
!     WATLEV: input  water levels                                         30.70
!
      REAL      BOTLEV(*), WATLEV(*)                                      30.70
!
!     FOUND : output  parameter indicating whether command
!                     being processed is found (value True)
!                     or not (False)
!
      LOGICAL   FOUND
!
!     Local variables

      REAL :: XPFR, YPFR, XLENFR, YLENFR                                  40.31

      TYPE(OPSDAT), POINTER :: OPSTMP, ROPS                               40.31

      TYPE XYPT                                                           40.31
        REAL                :: X, Y, XQ, YQ
        TYPE(XYPT), POINTER :: NEXTXY
      END TYPE XYPT

      TYPE(XYPT), TARGET  :: FRST                                         40.31
      TYPE(XYPT), POINTER :: CURR, TMP                                    40.31

      INTEGER, ALLOCATABLE :: VM(:)                                       40.80
      REAL, ALLOCATABLE :: XG(:), YG(:)                                   40.80
      CHARACTER (LEN=80) :: BASENM                                        40.80

!  8. Subroutines used
!
!     command reading routines
!     (all Ocean Pack)

      LOGICAL :: STPNOW                                                   34.01
      LOGICAL :: EQREAL ! if True the two (real) arguments are equal      33.09

!  9. Subroutines calling
!
!     SPROUT
!
! 10. Error messages
!
!     ---
!
! 13. Source text
!
      LOGICAL   PP                                                        30.72
      INTEGER   IERR, SIRAY                                               30.81 30.72
      CHARACTER PSNAME*16, STYPE*1, PRNAME*16                             40.31 30.21
      LOGICAL   KEYWIS, BOTDEP                                            30.70
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SWREPS')
!
!   --------------------------------------------------------------------
!   FRAME   'sname'  [xpfr] [ypfr] [alpfr] [xlenfr] [ylenfr]          &
!           [mxfr] [myfr]
!   --------------------------------------------------------------------
!
      IF (KEYWIS ('FRA')) THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (FRA) in combination with'//   32.02
     &                   ' 1D-computation')                               32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
!         ver 30.20: names of input variables changed, order of data changed
          ALLOCATE(OPSTMP)                                                40.31
          CALL INCSTR ('SNAME',PSNAME,'REQ',' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
          OPSTMP%PSNAME = PSNAME                                          40.31
          CALL READXY ('XPFR', 'YPFR', XPFR, YPFR, 'REQ', 0., 0.)         40.31 30.20
          OPSTMP%OPR(1) = XPFR                                            40.31
          OPSTMP%OPR(2) = YPFR                                            40.31
          CALL INREAL('ALPFR',ALPK,'REQ',0.)                              30.20
          IF (KSPHER.GT.0 .AND. .NOT.EQREAL(ALPK,0.)) CALL MSGERR (2,
     &          '[alpfr] must be 0 with spherical coordinates')           33.09
          CALL INREAL('XLENFR', XLENFR,'REQ',0.)                          40.31 30.20
          CALL INREAL('YLENFR', YLENFR,'REQ',0.)                          40.31 30.20
          OPSTMP%OPR(3) = XLENFR                                          40.31
          OPSTMP%OPR(4) = YLENFR                                          40.31
          OPSTMP%OPR(5) = PI2 * (ALPK/360.-NINT(ALPK/360.))
!           ***** the user gives number of meshes along each side *****
!           ***** program uses the number of points               *****
          CALL ININTG ('MXFR',MXK,'STA',20)                               30.20
          CALL ININTG ('MYFR',MYK,'STA',20)                               30.20
          OPSTMP%PSTYPE = 'F'                                             40.31
          OPSTMP%OPI(1) = MXK+1                                           40.31
          OPSTMP%OPI(2) = MYK+1                                           40.31
          ALLOCATE(OPSTMP%XP(0))                                          40.31
          ALLOCATE(OPSTMP%YP(0))                                          40.31
          NULLIFY(OPSTMP%NEXTOPS)                                         40.31
          IF ( .NOT.LOPS ) THEN                                           40.31
             FOPS = OPSTMP                                                40.31
             COPS => FOPS                                                 40.31
             LOPS = .TRUE.                                                40.31
          ELSE                                                            40.31
             COPS%NEXTOPS => OPSTMP                                       40.31
             COPS => OPSTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.02
      ENDIF
!
!   ------------------------------------------------------------------
!   GROUP   'sname'  SUBGRID [ix1] [ix2] [iy1] [iy2]
!   ------------------------------------------------------------------
!
      IF (KEYWIS('GROUP') .OR. KEYWIS ('SUBG')) THEN                      970221
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (GROUP) in combination'//      32.02
     &                   ' with 1D-computation')                          32.02
          GOTO 800                                                        32.02
        ELSEIF (OPTG.EQ.5) THEN                                           40.80
          CALL MSGERR(2,
     &              ' Keyword GROUP not supported in unstructured grid')  40.80
          GOTO 800                                                        40.80
        ELSE                                                              32.02
!         mod 970221: GROUP is introduced as a new command instead of
!         an option SUBG within the Frame command
          ALLOCATE(OPSTMP)                                                40.31
          CALL INCSTR ('SNAME',PSNAME,'REQ',' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
          OPSTMP%PSNAME = PSNAME                                          40.31
          CALL INKEYW ('STA', ' ')
          CALL IGNORE ('SUBG')                                            970221
          CALL ININTG ('IX1', IX1, 'REQ', 0)
          CALL ININTG ('IX2', IX2, 'REQ', 0)
          CALL ININTG ('IY1', IY1, 'REQ', 0)
          CALL ININTG ('IY2', IY2, 'REQ', 0)
          IF (IX1 .LT. 0 .OR. IX2 .GT. MXCGL-1 .OR. IX1 .GT. IX2 .OR.
     &        IY1 .LT. 0 .OR. IY2 .GT. MYCGL-1 .OR. IY1 .GT. IY2) THEN
            CALL MSGERR (3, 'Check corners of GROUP (SUBGRID) command')
            CALL MSGERR (3, ' .........the values should be.........')
            CALL MSGERR (3, 'ix1<ix2 and both between 0 and MXC')
            CALL MSGERR (3, 'iy1<iy2 and both between 0 and MYC')
          ENDIF
!
          IF (OPTG .EQ. 3) THEN
!             *** If the comput grid is curvilinear then the next    ***
!             *** quantities are stored : 'H' ,FLOAT(IX2), FLOAT(IY2)***
!             *** FLOAT(IX1) ,FLOAT(IY1) , 0 ,MXK+1 ,MYK+1           ***
!              *** Here frame type H is introduce, means that regular***
!              *** frame is required from a curvilinear compt. grid  ***
            OPSTMP%PSTYPE = 'H'                                           40.31
            MXK = IX2-IX1
            MYK = IY2-IY1
            OPSTMP%OPR(1) = FLOAT(IX2)                                    40.31
            OPSTMP%OPR(2) = FLOAT(IY2)                                    40.31
            OPSTMP%OPR(3) = FLOAT(IX1)                                    40.31
            OPSTMP%OPR(4) = FLOAT(IY1)                                    40.31
            OPSTMP%OPR(5) = 0.                                            40.31
            OPSTMP%OPI(1) = MXK+1                                         40.31
            OPSTMP%OPI(2) = MYK+1                                         40.31
          ELSE IF (OPTG .EQ. 1) THEN
            OPSTMP%PSTYPE = 'F'                                           40.31
            IF (IX1.NE.IX2) THEN
              OPSTMP%OPR(3) = (IX2-IX1)*DX                                40.31
            ELSE
              OPSTMP%OPR(3) = 0.01                                        40.31
            ENDIF
            IF (IY1.NE.IY2) THEN
              OPSTMP%OPR(4) = (IY2-IY1)*DY                                40.31
            ELSE
              OPSTMP%OPR(4) = 0.01                                        40.31
            ENDIF
            OPSTMP%OPR(1) = XPC + IX1*DX*COSPC - IY1*DY*SINPC             40.31
            OPSTMP%OPR(2) = YPC + IX1*DX*SINPC + IY1*DY*COSPC             40.31
            OPSTMP%OPR(5) = ALPC                                          40.31
            MXK = IX2-IX1
            MYK = IY2-IY1
            OPSTMP%OPI(1) = MXK+1                                         40.31
            OPSTMP%OPI(2) = MYK+1                                         40.31
            IF (ITEST .GE. 20 .OR. INTES .GE. 10)
     &        WRITE (PRINTF, 6020) (OPSTMP%OPR(II), II=2,6)               40.31
 6020       FORMAT (' Subgrid parms.', 6(1X,E12.4))
          ENDIF
          ALLOCATE(OPSTMP%XP(0))                                          40.31
          ALLOCATE(OPSTMP%YP(0))                                          40.31
          NULLIFY(OPSTMP%NEXTOPS)                                         40.31
          IF ( .NOT.LOPS ) THEN                                           40.31
             FOPS = OPSTMP                                                40.31
             COPS => FOPS                                                 40.31
             LOPS = .TRUE.                                                40.31
          ELSE                                                            40.31
             COPS%NEXTOPS => OPSTMP                                       40.31
             COPS => OPSTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.01
      ENDIF
!
!   ------------------------------------------------------------------
!   CURVE   'sname'  [xp1] [yp1]   < [int]  [xp]  [yp] >
!   ------------------------------------------------------------------
!

      IF (KEYWIS ('CURV')) THEN
        ALLOCATE(OPSTMP)                                                  40.31
        CALL INCSTR('SNAME',PSNAME,'REQ',' ')
        IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
        OPSTMP%PSNAME = PSNAME                                            40.31
        OPSTMP%PSTYPE = 'C'                                               40.31
        MIP  = 0
        OPSTMP%MIP = MIP                                                  40.31
!       ***** first point of a curve *****
   30   CALL NWLINE
        IF (STPNOW()) RETURN                                              34.01
        CALL READXY ('XP1', 'YP1', XP, YP, 'REQ', 0., 0.)
        FRST%X = XP                                                       40.31
        FRST%Y = YP                                                       40.31
        NULLIFY(FRST%NEXTXY)                                              40.31
        CURR => FRST                                                      40.31
        MIP = 1
!       ***** interval and next corner point *****
   33   CALL ININTG ('INT',INTV,'REP',-1)
        IF (INTV .NE. -1) THEN
          IF (INTV .LE. 0) THEN
             CALL MSGERR (2,'INT is negative or zero')
             INTV = 1
          ENDIF
          XP1 = XP
          YP1 = YP
          CALL READXY ('XP', 'YP', XP, YP, 'REQ', 0., 0.)                 40.03
          IF (ITEST .GE. 200 .OR. INTES .GE. 20) THEN                     30.21
            WRITE(PRINTF, 31) PSNAME
   31       FORMAT ('COORDINATES OF OUTPUT POINTS FOR CURVE  : ', A)
          ENDIF
          DO 36  JJ=1,INTV
            MIP = MIP+1
            ALLOCATE(TMP)                                                 40.31
            TMP%X = XP1+REAL(JJ)*(XP-XP1)/REAL(INTV)                      40.31
            TMP%Y = YP1+REAL(JJ)*(YP-YP1)/REAL(INTV)                      40.31
            NULLIFY(TMP%NEXTXY)                                           40.31
            CURR%NEXTXY => TMP                                            40.31
            CURR => TMP                                                   40.31
   36     CONTINUE
          GOTO 33
        ENDIF
        ALLOCATE(OPSTMP%XP(MIP))                                          40.31
        ALLOCATE(OPSTMP%YP(MIP))                                          40.31
        CURR => FRST                                                      40.31
        DO JJ = 1, MIP                                                    40.31
           OPSTMP%XP(JJ) = CURR%X                                         40.31
           OPSTMP%YP(JJ) = CURR%Y                                         40.31
           IF (ITEST .GE. 200 .OR. INTES .GE. 50) THEN                    40.31
              WRITE(PRINTF,32) JJ, CURR%X, CURR%Y                         40.31
   32         FORMAT(' POINT(',I4,')','  (IX,IY) -> ',2F10.2)             40.31
           ENDIF                                                          40.31
           CURR => CURR%NEXTXY                                            40.31
        END DO                                                            40.31
        DEALLOCATE(TMP)                                                   40.31
!       ***** store number of points of the curve *****
        OPSTMP%MIP = MIP                                                  40.31
        IF (MIP .EQ. 0) CALL MSGERR(1,'No output points found')
        NULLIFY(OPSTMP%NEXTOPS)                                           40.31
        IF ( .NOT.LOPS ) THEN                                             40.31
           FOPS = OPSTMP                                                  40.31
           COPS => FOPS                                                   40.31
           LOPS = .TRUE.                                                  40.31
        ELSE                                                              40.31
           COPS%NEXTOPS => OPSTMP                                         40.31
           COPS => OPSTMP                                                 40.31
        END IF                                                            40.31
        GOTO 800
      ENDIF
!
!   ------------------------------------------------------------------
!   POINTS  'sname'  < [xp]  [yp]  >     |    FILE 'fname'
!   ------------------------------------------------------------------
!
      IF (KEYWIS ('POIN')) THEN
        ALLOCATE(OPSTMP)                                                  40.31
        CALL INCSTR('SNAME',PSNAME,'REQ',' ')
        IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
        OPSTMP%PSNAME = PSNAME                                            40.31
        OPSTMP%PSTYPE = 'P'                                               40.31
        MIP  = 0
        OPSTMP%MIP = MIP                                                  40.31
        CALL INKEYW ('STA', ' ')
        IF (KEYWIS('FILE')) THEN
          IOSTAT = 0
          NDS    = 0
          PP     = .TRUE.
          CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
          CALL FOR (NDS, FILENM, 'OF', IOSTAT)                            10.31
          IF (STPNOW()) RETURN                                            34.01
        ELSE
          PP = .FALSE.
        ENDIF
        FRST%X = 0.                                                       40.31
        FRST%Y = 0.                                                       40.31
        NULLIFY(FRST%NEXTXY)                                              40.31
        CURR => FRST                                                      40.31
        DO
          IF (PP) THEN
            IERR = 0                                                      10.18
            CALL REFIXY (NDS, XP, YP, IERR)                               10.18
            IF (IERR.EQ.-1) GOTO 47
            IF (IERR.EQ.-2) THEN
              CALL MSGERR (2, 'Error reading point coord. from file')
              GOTO 800
            ENDIF
          ELSE
            CALL READXY ('XP', 'YP', XP, YP, 'REP', -1.E10, -1.E10)
            IF (XP .LT. -0.9E10) GOTO 47
          ENDIF
          MIP = MIP+1
          ALLOCATE(TMP)                                                   40.31
          TMP%X = XP                                                      40.31
          TMP%Y = YP                                                      40.31
          NULLIFY(TMP%NEXTXY)                                             40.31
          CURR%NEXTXY => TMP                                              40.31
          CURR => TMP                                                     40.31
        ENDDO
  47    ALLOCATE(OPSTMP%XP(MIP))                                          40.31
        ALLOCATE(OPSTMP%YP(MIP))                                          40.31
        CURR => FRST%NEXTXY                                               40.31
        DO JJ = 1, MIP                                                    40.31
           OPSTMP%XP(JJ) = CURR%X                                         40.31
           OPSTMP%YP(JJ) = CURR%Y                                         40.31
           CURR => CURR%NEXTXY                                            40.31
        END DO                                                            40.31
        DEALLOCATE(TMP)                                                   40.31
!       ***** store number of output points *****
        OPSTMP%MIP = MIP                                                  40.31
        IF (MIP .EQ. 0) CALL MSGERR (2, 'No output points found')         10.32
        NULLIFY(OPSTMP%NEXTOPS)                                           40.31
        IF ( .NOT.LOPS ) THEN                                             40.31
           FOPS = OPSTMP                                                  40.31
           COPS => FOPS                                                   40.31
           LOPS = .TRUE.                                                  40.31
        ELSE                                                              40.31
           COPS%NEXTOPS => OPSTMP                                         40.31
           COPS => OPSTMP                                                 40.31
        END IF                                                            40.31
        GOTO 800
      ENDIF
!
!   -------------------------------------------------------------------
!   RAY     'rname'  [xp1] [yp1] [xq1] [yq1]                       &
!          <  [int]  [xp]  [yp]  [xq]  [yq]  >
!   -------------------------------------------------------------------
!
      IF (KEYWIS ('RAY'))  THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (RAY) in combination'//        32.02
     &                   ' with 1D-computation')                          32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
          ALLOCATE(OPSTMP)                                                40.31
          CALL INCSTR('RNAME',PSNAME,'REQ',' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'RNAME is too long')
          OPSTMP%PSNAME = PSNAME                                          40.31
          OPSTMP%PSTYPE = 'R'                                             40.31
          MIP  = 1
          OPSTMP%MIP = MIP                                                40.31
!         first ray
          CALL NWLINE
          IF (STPNOW()) RETURN                                            34.01
          CALL READXY ('XP1', 'YP1', XP, YP, 'REQ', 0., 0.)
          CALL READXY ('XQ1', 'YQ1', XQ, YQ, 'REQ', 0., 0.)
          FRST%X  = XP                                                    40.31
          FRST%Y  = YP                                                    40.31
          FRST%XQ = XQ                                                    40.31
          FRST%YQ = YQ                                                    40.31
          NULLIFY(FRST%NEXTXY)                                            40.31
          CURR => FRST                                                    40.31
!         following rays
  110     CALL ININTG ('INT',INTD,'REP',-1)
          IF (INTD .NE. -1) THEN
            IF (INTD .LE. 0) THEN
              CALL MSGERR(2, 'INT negative or zero')
              INTD = 1
            ENDIF
            XP1 = XP
            YP1 = YP
            XQ1 = XQ
            YQ1 = YQ
            CALL READXY ('XP', 'YP', XP, YP, 'REQ', 0., 0.)
            CALL READXY ('XQ', 'YQ', XQ, YQ, 'REQ', 0., 0.)
            DO 115 JJ=1,INTD
              MIP = MIP+1
              ALLOCATE(TMP)                                               40.31
              TMP%X  = XP1 + REAL(JJ)*(XP-XP1)/REAL(INTD)                 40.31
              TMP%Y  = YP1 + REAL(JJ)*(YP-YP1)/REAL(INTD)                 40.31
              TMP%XQ = XQ1 + REAL(JJ)*(XQ-XQ1)/REAL(INTD)                 40.31
              TMP%YQ = YQ1 + REAL(JJ)*(YQ-YQ1)/REAL(INTD)                 40.31
              NULLIFY(TMP%NEXTXY)                                         40.31
              CURR%NEXTXY => TMP                                          40.31
              CURR => TMP                                                 40.31
  115       CONTINUE
            GOTO 110
          ENDIF
          ALLOCATE(OPSTMP%XP(MIP))                                        40.31
          ALLOCATE(OPSTMP%YP(MIP))                                        40.31
          ALLOCATE(OPSTMP%XQ(MIP))                                        40.31
          ALLOCATE(OPSTMP%YQ(MIP))                                        40.31
          CURR => FRST                                                    40.31
          DO JJ = 1, MIP                                                  40.31
             OPSTMP%XP(JJ) = CURR%X                                       40.31
             OPSTMP%YP(JJ) = CURR%Y                                       40.31
             OPSTMP%XQ(JJ) = CURR%XQ                                      40.31
             OPSTMP%YQ(JJ) = CURR%YQ                                      40.31
             CURR => CURR%NEXTXY                                          40.31
          END DO                                                          40.31
          DEALLOCATE(TMP)                                                 40.31
!
!         ***** termination *****
          OPSTMP%MIP = MIP                                                40.31
          IF (MIP .EQ. 1) CALL MSGERR (1,'Only one ray is defined')
          NULLIFY(OPSTMP%NEXTOPS)                                         40.31
          IF ( .NOT.LOPS ) THEN                                           40.31
             FOPS = OPSTMP                                                40.31
             COPS => FOPS                                                 40.31
             LOPS = .TRUE.                                                40.31
          ELSE                                                            40.31
             COPS%NEXTOPS => OPSTMP                                       40.31
             COPS => OPSTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.02
      ENDIF
!
!   -------------------------------------------------------------------
!   ISOline 'sname'  'rname'  DEPTH / BOTTOM [dep]
!   -------------------------------------------------------------------
!
      IF (KEYWIS ('ISO')) THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (ISO) in combination'//        32.02
     &                   ' with 1D-computation')                          32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
          ALLOCATE(OPSTMP)                                                40.31
          CALL INCSTR ('SNAME',PSNAME,'REQ',' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
          OPSTMP%PSNAME = PSNAME                                          40.31
          CALL INCSTR ('RNAME', PRNAME, 'REQ', ' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'RNAME is too long')
          CALL INKEYW ('STA', 'DEP')
          IF (KEYWIS ('BOT')) THEN
            BOTDEP = .TRUE.                                               30.70
          ELSE
            CALL IGNORE ('DEP')
            BOTDEP = .FALSE.                                              30.70
            IF (DYNDEP) CALL MSGERR (2,'depths will vary with time')      40.00
          ENDIF
          CALL INREAL ('DEP',DP,'REP',-1.E10)
          ROPS => FOPS                                                    40.31
          DO                                                              40.31
            IF (ROPS%PSNAME.EQ.PRNAME) EXIT                               40.31
            IF (.NOT.ASSOCIATED(ROPS%NEXTOPS)) THEN                       40.31
               CALL MSGERR(2,'Set of rays not defined')                   40.31
               GOTO 800                                                   40.31
            END IF                                                        40.31
            ROPS => ROPS%NEXTOPS                                          40.31
          END DO                                                          40.31
          STYPE = ROPS%PSTYPE                                             40.31
          IF (STYPE .NE. 'R') THEN                                        40.31
             CALL MSGERR                                                  40.31
     &          (2,'Ray name set assigned to set of output locations')    40.31
             GOTO 800                                                     40.31
          END IF                                                          40.31
          MIPR  = ROPS%MIP                                                40.31
          OPSTMP%PSTYPE = 'C'                                             40.31
          MIP  = 0
          OPSTMP%MIP = MIP                                                40.31
          FRST%X = 0.                                                     40.31
          FRST%Y = 0.                                                     40.31
          NULLIFY(FRST%NEXTXY)                                            40.31
          CURR => FRST                                                    40.31
          DO 125 IK=1,MIPR
             XP   = ROPS%XP(IK)                                           40.31
             YP   = ROPS%YP(IK)                                           40.31
             XQ   = ROPS%XQ(IK)                                           40.31
             YQ   = ROPS%YQ(IK)                                           40.31
             II   = SIRAY (DP, XP, YP, XQ, YQ, XX, YY, BOTDEP,            30.70
     &                     BOTLEV, WATLEV)                                30.70
             IF (II.EQ.0) THEN
                WRITE (PRINTF, 6120) DP, XP+XOFFS, YP+YOFFS,              30.70
     &                                   XQ+XOFFS, YQ+YOFFS               30.70
 6120           FORMAT(' No point with depth ',F5.2,
     &                 ' is found in ray :',4F10.2)
             ELSE
                MIP = MIP+1
                ALLOCATE(TMP)                                             40.31
                TMP%X = XX                                                40.31
                TMP%Y = YY                                                40.31
                NULLIFY(TMP%NEXTXY)                                       40.31
                CURR%NEXTXY => TMP                                        40.31
                CURR => TMP                                               40.31
             ENDIF
  125     CONTINUE
          ALLOCATE(OPSTMP%XP(MIP))                                        40.31
          ALLOCATE(OPSTMP%YP(MIP))                                        40.31
          CURR => FRST%NEXTXY                                             40.31
          DO IK = 1, MIP                                                  40.31
             OPSTMP%XP(IK) = CURR%X                                       40.31
             OPSTMP%YP(IK) = CURR%Y                                       40.31
             CURR => CURR%NEXTXY                                          40.31
          END DO                                                          40.31
          DEALLOCATE(TMP)                                                 40.31
          IF (MIP.EQ.0) CALL MSGERR
     &               (2, 'No points with valid depth found')
!             ***** store number of points of the curve *****
          OPSTMP%MIP = MIP                                                40.31
          NULLIFY(OPSTMP%NEXTOPS)                                         40.31
          IF ( .NOT.LOPS ) THEN                                           40.31
             FOPS = OPSTMP                                                40.31
             COPS => FOPS                                                 40.31
             LOPS = .TRUE.                                                40.31
          ELSE                                                            40.31
             COPS%NEXTOPS => OPSTMP                                       40.31
             COPS => OPSTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.02
      ENDIF
!
!   -------------------------------------------------------------------
!                    | [xpn] [ypn] [alpn] [xlenn] [ylenn] [mxn] [myn]
!   NGRID  'sname'  <
!                    | UNSTRUCtured / -> TRIAngle \
!                                   \    EASYmesh / 'fname'
!   -------------------------------------------------------------------
!
      IF (KEYWIS ('NGR')) THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (NGR) in combination'//        32.02
     &                   ' with 1D-computation')                          32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
!         ver 30.20: names changed, order changed
          ALLOCATE(OPSTMP)                                                40.31
          CALL INCSTR('SNAME',PSNAME,'REQ',' ')
          IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
          OPSTMP%PSNAME = PSNAME                                          40.31
          OPSTMP%PSTYPE = 'N'                                             40.31
          CALL INKEYW ('STA', ' ')
          IF (KEYWIS('UNSTRUC')) THEN                                     40.80
             PP = .TRUE.
             CALL INKEYW('STA','TRIA')
             IF (KEYWIS('EASY')) THEN
                IOSTAT = 0
                NDS    = 0
                CALL INCSTR ('FNAME', BASENM, 'REQ', ' ')
                FILENM = TRIM(BASENM)//'.n'
                CALL FOR (NDS, FILENM, 'OF', IOSTAT)
                IF (STPNOW()) RETURN
!
!               --- read first line to determine number of vertices
!
                READ(NDS, *, END=950, ERR=910) NVTX
                IF(.NOT.ALLOCATED(XG)) ALLOCATE(XG(NVTX), STAT = ISTAT)
                IF ( ISTAT == 0 ) THEN
                   IF(.NOT.ALLOCATED(YG)) ALLOCATE(YG(NVTX), STAT=ISTAT)
                ENDIF
                IF ( ISTAT == 0 ) THEN
                   IF(.NOT.ALLOCATED(VM)) ALLOCATE(VM(NVTX), STAT=ISTAT)
                ENDIF
                IF ( ISTAT /= 0 ) THEN
                   CALL MSGERR ( 4,
     &             'Allocation problem in SWREPS: array XG, YG or VM ' )
                   RETURN
                ENDIF
!
!               --- read coordinates of vertices and boundary marker
!
                DO KK = 1, NVTX
                   READ(NDS, 100, END=950, ERR=910) XG(I),YG(I),VM(I)
                ENDDO
!
!               --- close file <name>.n
!
                CLOSE(NDS)
!
             ELSE
                CALL IGNORE('TRIA')
                IOSTAT = 0
                NDS    = 0
                CALL INCSTR ('FNAME', BASENM, 'REQ', ' ')
                FILENM = TRIM(BASENM)//'.node'
                CALL FOR (NDS, FILENM, 'OF', IOSTAT)
                IF (STPNOW()) RETURN
!
!               --- read first line to determine number of vertices
!
                READ(NDS, *, END=950, ERR=910) NVTX, NDIM, NATTR, NBMARK
                IF(.NOT.ALLOCATED(XG)) ALLOCATE(XG(NVTX), STAT = ISTAT)
                IF ( ISTAT == 0 ) THEN
                   IF(.NOT.ALLOCATED(YG)) ALLOCATE(YG(NVTX), STAT=ISTAT)
                ENDIF
                IF ( ISTAT == 0 ) THEN
                   IF(.NOT.ALLOCATED(VM)) ALLOCATE(VM(NVTX), STAT=ISTAT)
                ENDIF
                IF ( ISTAT /= 0 ) THEN
                   CALL MSGERR ( 4,
     &             'Allocation problem in SWREPS: array XG, YG or VM ' )
                   RETURN
                ENDIF
!
!               --- check if boundary marker has been specified
!
                IF ( NBMARK == 0 ) THEN
                   CALL MSGERR ( 4,
     &         'boundary marker for vertices/faces must be specified ' )
                   RETURN
                ENDIF
!
!               --- read coordinates of vertices and boundary marker
!
                IF ( NATTR == 0 ) THEN
                   DO KK = 1, NVTX
                      READ(NDS, *, END=950, ERR=910) I,XG(I),YG(I),VM(I)
                   ENDDO
                ELSE
                   DO KK = 1, NVTX
                      READ(NDS, *, END=950, ERR=910) I,XG(I),YG(I),RDUM,
     &                                               VM(I)
                   ENDDO
                ENDIF
!
!               --- close file <name>.node
!
                CLOSE(NDS)
!
             ENDIF
          ELSE
             PP = .FALSE.
          ENDIF
          IF (.NOT.PP) THEN                                               40.80
!         structured grid
             CALL READXY ('XPN', 'YPN', XPCN, YPCN , 'REQ', 0., 0.)       30.20
             CALL INREAL('ALPN',ALPCN,'REQ',0.)                           30.20
             CALL INREAL('XLENN',XNLEN,'REQ',0.)                          30.20
             CALL INREAL('YLENN',YNLEN,'REQ',0.)                          30.20
             ALTNP = ALPCN / 360.
             ANG = PI2 * (ALTNP - NINT(ALTNP))
!            estimate step size for output                                40.00
             IF (OPTG.EQ.1) THEN                                          40.13
               COSA2 = (COS(ANG-ALPC))**2                                 40.00
               SINA2 = (SIN(ANG-ALPC))**2                                 40.00
               DXN = DX*COSA2 + DY*SINA2                                  40.00
               DYN = DX*SINA2 + DY*COSA2                                  40.00
             ELSEIF (OPTG.EQ.3) THEN                                      40.80 40.13
!              curvilinear grid, DXN and DYN are average step size        40.13
               DXN = (XCLEN+YCLEN)/REAL(MXCGL+MYCGL)                      40.31 40.13
               DYN = DXN                                                  40.13
             ELSEIF (OPTG.EQ.5) THEN                                      40.80
!              unstructured grid, DXN and DYN are average grid size       40.80
               DXN = 0.5*(mingsiz+maxgsiz)                                40.80
               DYN = DXN                                                  40.80
             ENDIF                                                        40.13
             CALL ININTG ('MXN',MXN,'STA',MAX(1,NINT(XNLEN/DXN)))         40.00
             CALL ININTG ('MYN',MYN,'STA',MAX(1,NINT(YNLEN/DYN)))         40.00
             MIP = 0
             ALLOCATE(OPSTMP%XP(2*(MXN+MYN)))                             40.31
             ALLOCATE(OPSTMP%YP(2*(MXN+MYN)))                             40.31
             XF=XPCN
             YF=YPCN
!            *****   start to calculate the positions       ********
!            *****  of the boundary point in the four sides ********
             DO 50 I=1,4
               INTE=MXN
               ALON=XNLEN
               ANGLE=ANG+PI2*(90.*REAL(I-1))/360.
               IF (I .EQ. 2 .OR. I .EQ. 4) THEN
                 INTE=MYN
                 ALON=YNLEN
               ENDIF
               XI=XF
               YI=YF
               XF=XI+ALON*COS(ANGLE)
               YF=YI+ALON*SIN(ANGLE)
               DO 55 KK=1,INTE
                 MIP=MIP+1
                 OPSTMP%XP(MIP)=XI+REAL(KK)*(XF-XI)/REAL(INTE)
                 OPSTMP%YP(MIP)=YI+REAL(KK)*(YF-YI)/REAL(INTE)
   55          CONTINUE
   50        CONTINUE
!                ***** store number of points *****
             OPSTMP%MIP = MIP                                             40.31
             OPSTMP%OPR(1) = XNLEN                                        40.31
             OPSTMP%OPR(2) = YNLEN                                        40.31
             OPSTMP%OPR(3) = XPCN                                         40.31
             OPSTMP%OPR(4) = YPCN                                         40.31
             OPSTMP%OPR(5) = PI2 * (ALPCN/360.-NINT(ALPCN/360.))          40.31
             OPSTMP%OPI(1) = MXN                                          40.31
             OPSTMP%OPI(2) = MYN                                          40.31
          ELSE                                                            40.80
!         unstructured grid
!
             MIP = COUNT(MASK=VM/=0)
             ALLOCATE(OPSTMP%XP(MIP))
             ALLOCATE(OPSTMP%YP(MIP))
             I = 0
             DO KK = 1, NVTX
                IF ( VM(KK) /= 0 ) THEN
                   I = I + 1
                   OPSTMP%XP(I)= XG(KK)-XOFFS
                   OPSTMP%YP(I)= YG(KK)-YOFFS
                ENDIF
             ENDDO
             IF (ALLOCATED(XG)) DEALLOCATE(XG)
             IF (ALLOCATED(YG)) DEALLOCATE(YG)
             IF (ALLOCATED(VM)) DEALLOCATE(VM)
             OPSTMP%MIP = MIP
!            next variable will be used for checking grid in
!            present computational and/or bottom grid
             OPSTMP%OPR(1) = -999.
          ENDIF                                                           40.80
          IF (MIP .EQ. 0) CALL MSGERR(1,'No output points found')
          NULLIFY(OPSTMP%NEXTOPS)                                         40.31
          IF ( .NOT.LOPS ) THEN                                           40.31
             FOPS = OPSTMP                                                40.31
             COPS => FOPS                                                 40.31
             LOPS = .TRUE.                                                40.31
          ELSE                                                            40.31
             COPS%NEXTOPS => OPSTMP                                       40.31
             COPS => OPSTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
!
 910      INQUIRE (UNIT=NDS, NAME=FILENM)                                 40.80
          CALL MSGERR (4, 'error reading data from file '//FILENM )       40.80
          RETURN                                                          40.80
 950      INQUIRE (UNIT=NDS, NAME=FILENM)                                 40.80
          CALL MSGERR (4, 'unexpected end of file in file '//FILENM )     40.80
          RETURN                                                          40.80
!
        ENDIF                                                             32.02
      ENDIF
!     ---------------------------------------------------------
!     command not found:
      RETURN
 800  FOUND = .TRUE.
      RETURN
 100  FORMAT((6X,2E22.15,I3))                                             40.80
!*    end of subroutine SWREPS  **
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWREOQ ( FOUND )                                         40.31
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
      USE OUTP_DATA                                                       40.13
      USE M_PARALL                                                        40.31
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
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde (1D version)
!     34.01: Jeroen Adema
!     40.03, 40.13: Nico Booij
!     40.14: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.50         : option COORD added in command PLOT
!                     option STAR  added in command PLOT
!     32.02, Feb. 98: 1D version introduced
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Apr. 98: removed reference to commons KAART and KAR
!     30.81, Nov. 98: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     30.81, Jan. 99: Replaced variable TO by TO_ (because TO is a reserved
!                     word)
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Nov. 99: in case SPEC2D the value of MXOUTAR is increased by
!                     6*MIP
!     40.03, Mar. 00: NQUA increased in case of Isoline plot
!            Sep. 00: inconsistency with manual corrected
!     40.13, Mar. 01: option BLOCKed added in plot of problem points
!            Aug. 01: array for NESTOUT request extended to 20 (in view of SETUP)
!     40.13, Oct. 01: filenames are stored in array OUTP_FILES
!                     not any more in array containing output request parameters
!     40.14, Dec. 01: format for setup corrected.
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing POOL construction and HPGL funcationality
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Reading and processing of the output requests
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     FOUND : output parameter indicating whether command
!                    being processed is found (value True)
!                    or not (False)
!
      LOGICAL FOUND                                                       30.72
!
!  8. Subroutines used
!
!     command reading routines
!     (all Ocean Pack)
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SPROUT
!
! 10. Error messages
!
!       ---
!
! 13. Source text
!
      INTEGER IERR                                                        40.31
      CHARACTER PSNAME *16, STYPE *1, RTYPE *4                            40.31
      LOGICAL   KEYWIS                                                    40.31
      TYPE(ORQDAT), POINTER :: ORQTMP                                     40.31
      TYPE(ORQDAT), SAVE, POINTER :: CORQ                                 40.31
      LOGICAL, SAVE :: LORQ = .FALSE.                                     40.31
      TYPE AUXT                                                           40.31
        INTEGER             :: I
        REAL                :: R
        TYPE(AUXT), POINTER :: NEXTI
      END TYPE AUXT
      TYPE(AUXT), TARGET  :: FRST                                         40.31
      TYPE(AUXT), POINTER :: CURR, TMP                                    40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SWREOQ')
!
!   --------------------------------------------------------------------------
!   BLOCK   'sname'  HEADER / NOHEADER  'fname' (LAY-OUT [idla])             &
!          <  DSPR/HSIGN/DIR/PDIR/TDIR/TM01/RTM01/RTP/TM02/FSPR/DEPTH/VEL/   &
!             FRCOEFF/WIND/DISSIP/QB/TRANSP/FORCE/UBOT/URMS/WLEN/STEEPNESS/  &
!             DHSIGN/DRTM01/LEAK/TSEC/XP/YP/DIST/SETUP/TMM10/RTMM10/         &
!             TMBOT/QP/BFI/WATLEV/BOTLEV/TPS/DISBOT/DISSURF/DISWCAP/         &
!             GENE/GENW/REDI/REDQ/REDT/PROPA/PROPX/PROPT/PROPS/RADS|LWAVP >  &
!             ([unit]) (OUTPUT [tbegblk] [deltblk] SEC/MIN/HR/DAY)
!   --------------------------------------------------------------------------
!   BLO   block type output
      IF (KEYWIS ('BLO')) THEN
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (BLO) in combination'//        32.02
     &                   ' with 1D-computation')                          32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
          CALL SWNMPS (PSNAME, STYPE, MIP, IERR)                          40.31
          IF (IERR.NE.0) GOTO 800
          IF (STYPE.NE.'F' .AND. STYPE.NE.'H' .AND. STYPE.NE.'U' ) THEN   40.80 30.21
            CALL MSGERR(2,'Set of output locations is not correct type')
            GOTO 800
          ENDIF
!
!         output frame exists
!
          ALLOCATE(ORQTMP)                                                40.31
          NREOQ = NREOQ + 1                                               40.31
          IF (NREOQ.GT.MAX_OUTP_REQ) CALL MSGERR (2,                      40.31 40.13
     &    'too many output requests')                                     40.13

          IDLAO = 1
!
          CALL INKEYW ('REQ',' ')
          IF (KEYWIS('NOHEAD') .OR. KEYWIS ('FIL')) THEN                  30.20
            CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
            DFAC = 1.
            CALL INKEYW ('STA', ' ')
            IF (KEYWIS('LONG')) THEN
!             option disabled                                             40.13
              CALL MSGERR (2, 'option LONG disabled; use OUTP OPT')       40.13
            ENDIF                                                         40.13
            RTYPE = 'BLKD'
          ELSE
            CALL IGNORE ('HEAD')                                          30.20
            CALL IGNORE ('PAP')
            DFAC = -1.
            CALL INCSTR ('FNAME', FILENM, 'STA', ' ')                     30.20
            RTYPE = 'BLKP'
            IF ( INDEX( FILENM, '.MAT' ).NE.0 .OR.                        40.41 40.30
     &           INDEX (FILENM, '.mat' ).NE.0 ) THEN                      40.41 40.30
               CALL MSGERR(4,'No header allowed for Matlab files')        40.30
               RETURN                                                     40.30
            END IF                                                        40.30
          END IF
          IF (FILENM .EQ. ' ') THEN                                       24/FEB
            NREF = PRINTF
          ELSE
            NREF = 0
!           --- append node number to FILENM in case of                   40.30
!               parallel computing                                        40.30
            IF ( PARLL ) THEN                                             40.30
               ILPOS = INDEX ( FILENM, ' ' )-1                            40.30
               WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                    40.30
  33           FORMAT('-',I3.3)                                           40.30
            END IF                                                        40.30
!PUN            FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                 40.95
          ENDIF
          CALL INKEYW ('STA', ' ')
          IF (KEYWIS('LAY')) THEN
            CALL ININTG ('IDLA', IDLAO, 'REQ', 0)
            CALL INKEYW ('REQ', ' ')
            IF (IDLAO.NE.1 .AND. IDLAO.NE.3 .AND. IDLAO.NE.4)
     &        CALL MSGERR (2, 'Illegal value for IDLA')
          ENDIF
          ORQTMP%OQR(1) = -1.                                             40.31
          ORQTMP%OQR(2) = -1.                                             40.31
          ORQTMP%RQTYPE = RTYPE                                           40.31
          ORQTMP%PSNAME = PSNAME                                          40.31
          ORQTMP%OQI(1) = NREF                                            40.31
          ORQTMP%OQI(2) = NREOQ                                           40.31
          NVAR = 0                                                        40.31
          ORQTMP%OQI(3) = NVAR                                            40.31
          ORQTMP%OQI(4) = IDLAO                                           40.31
          OUTP_FILES(NREOQ) = FILENM                                      40.13
!
!         read types of output quantities
!
          FRST%I = 0                                                      40.31
          FRST%R = 0.                                                     40.31
          NULLIFY(FRST%NEXTI)                                             40.31
          CURR => FRST                                                    40.31
  70      CALL SVARTP (IVTYPE)
          IF (IVTYPE .EQ. 98) GOTO 91                                     30.00
          IF (IVTYPE .NE. 99) THEN
             CALL INREAL ('UNIT', DFAC, 'STA', -1.)
             IF (OVSVTY(IVTYPE).EQ.5) THEN
               CALL MSGERR (2,
     &         'Type of output not allowed for this quantity')
               WRITE (PRINTF, *) ' -> ', OVSNAM(IVTYPE)
             ELSE IF (IVTYPE .GT. 0) THEN
                NVAR = NVAR+1
                ALLOCATE(TMP)                                             40.31
                TMP%I = IVTYPE                                            40.31
                TMP%R = DFAC                                              40.31
                NULLIFY(TMP%NEXTI)                                        40.31
                CURR%NEXTI => TMP                                         40.31
                CURR => TMP                                               40.31
                IF (IVTYPE.EQ.6) IUBOTR = 1
                IF (IVTYPE.EQ.50 .AND. JPBOT.LE.1) THEN                   40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JPBOT  = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.52 .AND. JBOTLV.LE.1) THEN                  40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JBOTLV = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.54 .AND. JDSXB.LE.1) THEN                   40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JDSXB  = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.55 .AND. JDSXS.LE.1) THEN                   40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JDSXS  = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.56 .AND. JDSXW.LE.1) THEN                   40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JDSXW  = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.57 .AND. JDSXV.LE.1) THEN                   40.65
                   MCMVAR = MCMVAR+1                                      40.65
                   JDSXV  = MCMVAR                                        40.65
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.65
                IF (IVTYPE.EQ.60 .AND. JGENR.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JGENR  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.61 .AND. JGSXW.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JGSXW  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.62 .AND. JREDS.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JREDS  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.63 .AND. JRSXQ.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JRSXQ  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.64 .AND. JRSXT.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JRSXT  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.65 .AND. JTRAN.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JTRAN  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.66 .AND. JTSXG.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JTSXG  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.67 .AND. JTSXT.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JTSXT  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.68 .AND. JTSXS.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JTSXS  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.69 .AND. JRADS.LE.1) THEN                   40.85
                   MCMVAR = MCMVAR+1                                      40.85
                   JRADS  = MCMVAR                                        40.85
                   ALOCMP = .TRUE.                                        40.97
                ENDIF                                                     40.85
                IF (IVTYPE.EQ.7  .OR. IVTYPE.EQ.9  .OR.                   40.85
     &              IVTYPE.EQ.54 .OR. IVTYPE.EQ.55 .OR.                   40.85
     &              IVTYPE.EQ.56 .OR. IVTYPE.GE.60 ) LADDS = .TRUE.       40.85
             ENDIF
             GOTO 70
          ENDIF
 91       IF (NVAR.GT.0) THEN                                             40.31
             ALLOCATE(ORQTMP%IVTYP(NVAR))                                 40.31
             ALLOCATE(ORQTMP%FAC(NVAR))                                   40.31
             CURR => FRST%NEXTI                                           40.31
             DO JJ = 1, NVAR                                              40.31
                ORQTMP%IVTYP(JJ) = CURR%I                                 40.31
                ORQTMP%FAC  (JJ) = CURR%R                                 40.31
                CURR => CURR%NEXTI                                        40.31
             END DO                                                       40.31
             DEALLOCATE(TMP)                                              40.31
          END IF                                                          40.31
!
          IF (IVTYPE .EQ. 98) THEN                                        30.00
            IF (NSTATM.EQ.0) CALL MSGERR (3,
     &      'time information not allowed in stationary mode')
            NSTATM = 1
            CALL INCTIM (ITMOPT, 'TBEG', ORQTMP%OQR(1), 'REQ', 0.)        40.31 30.00
            CALL ININTV ('DELT', ORQTMP%OQR(2), 'REQ', 0.)                40.31 30.00
          ENDIF
!
          ORQTMP%OQI(3) = NVAR                                            40.31
          NULLIFY(ORQTMP%NEXTORQ)                                         40.31
          IF ( .NOT.LORQ ) THEN                                           40.31
             FORQ = ORQTMP                                                40.31
             CORQ => FORQ                                                 40.31
             LORQ = .TRUE.                                                40.31
          ELSE                                                            40.31
             CORQ%NEXTORQ => ORQTMP                                       40.31
             CORQ => ORQTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.02
      ENDIF
!   --------------------------------------------------------------------------
!   TABLE   'sname'  HEADER / NOHEADER / INDEXED 'fname'                     &
!          <  DSPR/HSIGN/DIR/PDIR/TDIR/TM01/RTM01/RTP/TM02/FSPR/DEPTH/VEL/   &
!             FRCOEFF/WIND/DISSIP/QB/TRANSP/FORCE/UBOT/URMS/WLEN/STEEPNESS/  &
!             DHSIGN/DRTM01/LEAK/TIME/TSEC/XP/YP/DIST/SETUP/TMM10/RTMM10/    &
!             TMBOT/QP/BFI/WATLEV/BOTLEV/TPS/DISBOT/DISSURF/DISWCAP/         &
!             GENE/GENW/REDI/REDQ/REDT/PROPA/PROPX/PROPT/PROPS/RADS|LWAVP >  &
!             ([unit]) (OUTPUT [tbegtbl] [delttbl] SEC/MIN/HR/DAY)
!   --------------------------------------------------------------------------
!   TABLE   output in the form of a table

      IF (KEYWIS ('TAB')) THEN
        CALL SWNMPS (PSNAME, STYPE, MIP, IERR)                            40.31
        IF (IERR.NE.0) GOTO 800
!
!       output points exist
!
        ALLOCATE(ORQTMP)                                                  40.31
        NREOQ = NREOQ + 1                                                 40.31
        IF (NREOQ.GT.MAX_OUTP_REQ) CALL MSGERR (2,                        40.31 40.13
     &    'too many output requests')                                     40.13
!
        CALL INKEYW ('STA','HEAD')                                        20.67
        IF (KEYWIS('NOHEAD') .OR. KEYWIS ('FIL')) THEN                    20.67
          RTYPE = 'TABD'
        ELSE IF (KEYWIS ('IND')) THEN                                     30.50
          RTYPE = 'TABI'
        ELSE IF (KEYWIS ('SWAN')) THEN                                    40.00
          RTYPE = 'TABS'
        ELSE IF (KEYWIS ('STAB')) THEN                                    40.00
          RTYPE = 'TABT'
        ELSE
          CALL IGNORE ('HEAD')                                            20.67
          CALL IGNORE ('PAP')
          RTYPE = 'TABP'
        END IF
        ORQTMP%OQR(1) = -1.                                               40.31
        ORQTMP%OQR(2) = -1.                                               40.31
        ORQTMP%RQTYPE = RTYPE                                             40.31
!       unit reference number NREF is 0, will be determined in output module
        CALL INCSTR ('FNAME', FILENM, 'STA', ' ')
        IF (FILENM .NE. '    ') THEN
          NREF = 0
!         --- append node number to FILENM in case of                     40.30
!             parallel computing                                          40.30
          IF ( PARLL ) THEN                                               40.30
             ILPOS = INDEX ( FILENM, ' ' )-1                              40.30
             WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                      40.30
          END IF                                                          40.30
!PUN          FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                   40.95
        ELSE
          NREF = PRINTF
        ENDIF
        ORQTMP%PSNAME = PSNAME                                            40.31
        ORQTMP%OQI(1) = NREF                                              40.31
        ORQTMP%OQI(2) = NREOQ                                             40.31
        OUTP_FILES(NREOQ) = FILENM                                        40.31 40.13
!
        NVAR = 0
        ORQTMP%OQI(3) = NVAR                                              40.31
!       read types of variables to be printed in the table
        FRST%I = 0                                                        40.31
        NULLIFY(FRST%NEXTI)                                               40.31
        CURR => FRST                                                      40.31
   80   CALL SVARTP (IVTYPE)
        IF (IVTYPE .EQ. 98) GOTO 90                                       30.00
        IF (IVTYPE .NE. 99) THEN
          IF (OVSVTY(IVTYPE).EQ.5) THEN
            CALL MSGERR (2,
     &      'Type of output not allowed for this quantity')
            WRITE (PRINTF, *) ' -> ', OVSNAM(IVTYPE)
          ELSE IF (IVTYPE .GT. 0) THEN
            NVAR = NVAR+1
            ALLOCATE(TMP)                                                 40.31
            TMP%I = IVTYPE                                                40.31
            NULLIFY(TMP%NEXTI)                                            40.31
            CURR%NEXTI => TMP                                             40.31
            CURR => TMP                                                   40.31
            IF (IVTYPE.EQ.18) IUBOTR = 1
            IF (IVTYPE.EQ.50 .AND. JPBOT.LE.1) THEN                       40.65
               MCMVAR = MCMVAR+1                                          40.65
               JPBOT  = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.52 .AND. JBOTLV.LE.1) THEN                      40.65
               MCMVAR = MCMVAR+1                                          40.65
               JBOTLV = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.54 .AND. JDSXB.LE.1) THEN                       40.65
               MCMVAR = MCMVAR+1                                          40.65
               JDSXB  = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.55 .AND. JDSXS.LE.1) THEN                       40.65
               MCMVAR = MCMVAR+1                                          40.65
               JDSXS  = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.56 .AND. JDSXW.LE.1) THEN                       40.65
               MCMVAR = MCMVAR+1                                          40.65
               JDSXW  = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.57 .AND. JDSXV.LE.1) THEN                       40.65
               MCMVAR = MCMVAR+1                                          40.65
               JDSXV  = MCMVAR                                            40.65
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.65
            IF (IVTYPE.EQ.60 .AND. JGENR.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JGENR  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.61 .AND. JGSXW.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JGSXW  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.62 .AND. JREDS.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JREDS  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.63 .AND. JRSXQ.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JRSXQ  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.64 .AND. JRSXT.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JRSXT  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.65 .AND. JTRAN.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JTRAN  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.66 .AND. JTSXG.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JTSXG  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.67 .AND. JTSXT.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JTSXT  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.68 .AND. JTSXS.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JTSXS  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.69 .AND. JRADS.LE.1) THEN                       40.85
               MCMVAR = MCMVAR+1                                          40.85
               JRADS  = MCMVAR                                            40.85
               ALOCMP = .TRUE.                                            40.97
            ENDIF                                                         40.85
            IF (IVTYPE.EQ.7  .OR. IVTYPE.EQ.9  .OR.                       40.85
     &          IVTYPE.EQ.54 .OR. IVTYPE.EQ.55 .OR.                       40.85
     &          IVTYPE.EQ.56 .OR. IVTYPE.GE.60 ) LADDS = .TRUE.           40.85
            CALL INKEYW ('STA', ' ')                                      40.00
            IF (KEYWIS('UNIT')) THEN
              CALL MSGERR (1, 'UNIT is ignored in this version')          40.00
            ENDIF
          ENDIF
          GOTO 80
        ENDIF
 90     IF (NVAR.GT.0) THEN                                               40.31
           ALLOCATE(ORQTMP%IVTYP(NVAR))                                   40.31
           CURR => FRST%NEXTI                                             40.31
           DO JJ = 1, NVAR                                                40.31
              ORQTMP%IVTYP(JJ) = CURR%I                                   40.31
              CURR => CURR%NEXTI                                          40.31
           END DO                                                         40.31
           DEALLOCATE(TMP)                                                40.31
        END IF                                                            40.31
        ALLOCATE(ORQTMP%FAC(0))                                           40.31

        IF (IVTYPE .EQ. 98) THEN                                          30.00
          IF (NSTATM.EQ.0) CALL MSGERR (3,
     &      'time information not allowed in stationary mode')
          NSTATM = 1
          CALL INCTIM (ITMOPT, 'TBEG', ORQTMP%OQR(1), 'REQ', 0.)          40.31 30.00
          CALL ININTV ('DELT', ORQTMP%OQR(2), 'REQ', 0.)                  40.31 30.00
        ENDIF
        ORQTMP%OQI(3) = NVAR                                              40.31
        NULLIFY(ORQTMP%NEXTORQ)                                           40.31
        IF ( .NOT.LORQ ) THEN                                             40.31
           FORQ = ORQTMP                                                  40.31
           CORQ => FORQ                                                   40.31
           LORQ = .TRUE.                                                  40.31
        ELSE                                                              40.31
           CORQ%NEXTORQ => ORQTMP                                         40.31
           CORQ => ORQTMP                                                 40.31
        END IF                                                            40.31
        GOTO 800
      ENDIF
!
!   PLOT    plot iso lines and/or vector fields
!
!   --------------------------------------------------------------------------

      IF (KEYWIS ('PLO')) THEN
        CALL MSGERR(2,'Keyword PLO... is no longer maintained')           40.31
        GOTO 800
      ENDIF
!
!   --------------------------------------------------------------------------
!   SPECout 'sname'  SPEC1D/SPEC2D  ABS/REL   'fname'                       &
!             (OUTPUT [tbegspc] [deltspc] SEC/MIN/HR/DAY)
!   --------------------------------------------------------------------------
!   SPEC   output of spectra

      IF (KEYWIS ('SPEC')) THEN
        CALL SWNMPS (PSNAME, STYPE, MIP, IERR)                            40.31
        IF (IERR.NE.0) GOTO 800
!
!       output points exist
!
        ALLOCATE(ORQTMP)                                                  40.31
        NREOQ = NREOQ + 1                                                 40.31
        IF (NREOQ.GT.MAX_OUTP_REQ) CALL MSGERR (2,                        40.31 40.13
     &    'too many output requests')                                     40.13
!
        CALL INKEYW ('STA', 'SPEC2D')                                     20.67
        IF (KEYWIS ('FS1D') .OR. KEYWIS('SPEC1D')) THEN                   20.67
          RTYPE  = 'SPE1'                                                 20.28
        ELSE
          CALL IGNORE ('SFD')                                             20.28
          CALL IGNORE ('SPEC2D')                                          20.67
          RTYPE  = 'SPEC'
        ENDIF
        CALL INKEYW ('STA', 'ABS')                                        40.03
        IF (KEYWIS ('REL')) THEN                                          40.03
          RTYPE(3:3)  = 'R'                                               20.28
        ELSE
          CALL IGNORE ('ABS')                                             40.03
        ENDIF
        ORQTMP%OQR(1) = -1.                                               40.31
        ORQTMP%OQR(2) = -1.                                               40.31
        ORQTMP%RQTYPE = RTYPE                                             40.31
        CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
        NREF = 0
!       --- append node number to FILENM in case of                       40.30
!           parallel computing                                            40.30
        IF ( PARLL ) THEN                                                 40.30
           ILPOS = INDEX ( FILENM, ' ' )-1                                40.30
           WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                        40.30
        END IF                                                            40.30
!PUN        FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                     40.95
        ORQTMP%PSNAME = PSNAME                                            40.31
        ORQTMP%OQI(1) = NREF                                              40.31
        ORQTMP%OQI(2) = NREOQ                                             40.31
        OUTP_FILES(NREOQ) = FILENM                                        40.31 40.13
!
        NVAR = 0                                                          40.31
        ORQTMP%OQI(3) = NVAR                                              40.31
        ALLOCATE(ORQTMP%IVTYP(0))                                         40.31
        ALLOCATE(ORQTMP%FAC(0))                                           40.31
!       read types of variables to be printed in the table
        CALL INKEYW ('STA', ' ')                                          30.00
        IF (KEYWIS ('OUT')) THEN                                          40.03
          IF (NSTATM.EQ.0) CALL MSGERR (3,
     &      'time information not allowed in stationary mode')
          NSTATM = 1
          CALL INCTIM (ITMOPT, 'TBEG', ORQTMP%OQR(1), 'REQ', 0.)          40.31 30.00
          CALL ININTV ('DELT', ORQTMP%OQR(2), 'REQ', 0.)                  40.31 30.00
          IF (NSTATM.EQ.0) CALL MSGERR (2,
     &                  'time input not allowed in stationary mode')
        ENDIF
        NULLIFY(ORQTMP%NEXTORQ)                                           40.31
        IF ( .NOT.LORQ ) THEN                                             40.31
           FORQ = ORQTMP                                                  40.31
           CORQ => FORQ                                                   40.31
           LORQ = .TRUE.                                                  40.31
        ELSE                                                              40.31
           CORQ%NEXTORQ => ORQTMP                                         40.31
           CORQ => ORQTMP                                                 40.31
        END IF                                                            40.31
        GOTO 800
      END IF

!   --------------------------------------------------------------------------
!   NESTout 'sname'  'fname'                                                &
!             (OUTPUT [tbegnst] [deltnst] SEC/MIN/HR/DAY)
!   --------------------------------------------------------------------------
!   NEST   output for nesting of models                        VER.       20.63

      IF (KEYWIS ('NEST')) THEN                                           40.00
!
!      ======================================================================
!
!       NESTout  'sname'  'fname'  &
!
!                                              | -> Sec  |
!                OUTput  [tbegnst]  [deltnst] <     MIn   >
!                                              |    HR   |
!                                              |    DAy  |
!
!      =======================================================================
!
        IF (ONED) THEN                                                    32.02
          CALL MSGERR (2,' Illegal keyword (NEST) in'//                   32.02
     &                   ' combination with 1D-computation')              32.02
          GOTO 800                                                        32.02
        ELSE                                                              32.02
          CALL SWNMPS (PSNAME, STYPE, MIP, IERR)                          40.31
          IF (IERR.NE.0) GOTO 800
          IF (STYPE .NE. 'N') THEN
            CALL MSGERR(2,'Set of output locations is not correct type')
            GOTO 800
          ENDIF
!         output points exist
          ALLOCATE(ORQTMP)                                                40.31
          NREOQ = NREOQ + 1                                               40.31
          IF (NREOQ.GT.MAX_OUTP_REQ) CALL MSGERR (2,                      40.31 40.13
     &    'too many output requests')                                     40.13
          ORQTMP%OQR(1) = -1.                                             40.31
          ORQTMP%OQR(2) = -1.                                             40.31
          RTYPE  = 'SPRC'                                                 40.00
          ORQTMP%RQTYPE = RTYPE                                           40.31
          CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
          NREF = 0
!         --- append node number to FILENM in case of                     40.30
!             parallel computing                                          40.30
          IF ( PARLL ) THEN                                               40.30
             ILPOS = INDEX ( FILENM, ' ' )-1                              40.30
             WRITE(FILENM(ILPOS+1:ILPOS+4),33) INODE                      40.30
          END IF                                                          40.30
!PUN          FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                   40.95
          ORQTMP%PSNAME = PSNAME                                          40.31
          ORQTMP%OQI(1) = NREF                                            40.31
          ORQTMP%OQI(2) = NREOQ                                           40.31
          OUTP_FILES(NREOQ) = FILENM                                      40.31 40.13
          NVAR = 0                                                        40.31
          ORQTMP%OQI(3) = NVAR                                            40.31
          ALLOCATE(ORQTMP%IVTYP(0))                                       40.31
          ALLOCATE(ORQTMP%FAC(0))                                         40.31
!
          CALL INKEYW ('STA', ' ')                                        30.00
          IF (KEYWIS ('OUT')) THEN                                        40.03
            IF (NSTATM.EQ.0) CALL MSGERR (3,
     &      'time information not allowed in stationary mode')
            NSTATM = 1
            CALL INCTIM (ITMOPT, 'TBEG', ORQTMP%OQR(1), 'REQ', 0.)        40.31 30.00
            CALL ININTV ('DELT', ORQTMP%OQR(2), 'REQ', 0.)                40.31 30.00
            IF (NSTATM.EQ.0) CALL MSGERR (2,
     &                  'time input not allowed in stationary mode')
          ENDIF
!
          NULLIFY(ORQTMP%NEXTORQ)                                         40.31
          IF ( .NOT.LORQ ) THEN                                           40.31
             FORQ = ORQTMP                                                40.31
             CORQ => FORQ                                                 40.31
             LORQ = .TRUE.                                                40.31
          ELSE                                                            40.31
             CORQ%NEXTORQ => ORQTMP                                       40.31
             CORQ => ORQTMP                                               40.31
          END IF                                                          40.31
          GOTO 800
        ENDIF                                                             32.02
      ENDIF
!     -------------------------------------------------------
!     command not found:
      RETURN
 800  FOUND = .TRUE.
      RETURN
!*    end of subroutine SWREOQ  **
      END
!***********************************************************************
!                                                                      *
      INTEGER FUNCTION SIRAY (DP, XP1, YP1, XP2, YP2, XX, YY, BOTDEP,     30.70
     &                        BOTLEV, WATLEV)                             30.70
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
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
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     00.00, Mar. 87: heading added, name of routine changed from
!                     IRAAI in SIRAY
!     30.72, Oct. 97: logical function EQREAL introduced for floating point
!                     comparisons
!     30.70, Nov. 97: changed into INTEGER function
!                     test output added
!                     arguments BOTDEP, BOTLEV, WATLEV added
!     40.03, Nov. 99: X2= etc. moved out of IF-ENDIF group
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Searching the first point on a ray where the depth is DP
!
!  3. METHOD
!
!     ---
!
!  4. PARAMETERLIST
!
!     DP      REAL   input    depth
!     XP1     REAL   input    X-coordinate start point of ray
!     YP1     REAL   input    Y-coordinate start point of ray
!     XP2     REAL   input    X-coordinate end point of ray
!     YP2     REAL   input    Y-coordinate end point of ray
!     XX      REAL   input    X-coordinate point with depth DP
!     YY      REAL   input    Y-coordinate point with depth DP
!
!  5. SUBROUTINES CALLING
!
!     SWREPS (SWAN/READ)
!
!  6. SUBROUTINES USED
!
!     SVALQI (SWAN/READ)
!
!  7. ERROR MESSAGES
!
!     ---
!
!  8. REMARKS
!
!     ---
!
!  9. STRUCTURE
!
!     ----------------------------------------------------------------
!     Give SIRAY initial value 0
!     Compute stepsize, raylength and number of steps along the ray
!     Compute bottom coordinates  of startpoint as number of meshes
!     Call SVALQI to interpolate depth in startpoint of ray
!     For every step along the ray do
!         Compute coordinates of the intermediate point in problem
!           grid and bottom grid
!         Call SVALQI to interpolate the depth for this point
!         If the required depth is in the interval, then
!             Compute coordinates of the point with depth DP
!             Set SIRAY 1
!         Else
!             Coordinates and depth at start of new interval are values
!               at end of old interval
!     ----------------------------------------------------------------
!  10. SOURCE TEXT
!
      LOGICAL   EQREAL, BOTDEP                                            30.72
      REAL      BOTLEV(*), WATLEV(*)                                      30.70
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SIRAY')
!
      SIRAY   = 0
      DIFDEP = 1e+10
      DSTEP  = MIN(DXG(1), DYG(1))
      RAYLEN = SQRT ((XP2-XP1)*(XP2-XP1) + (YP2-YP1)*(YP2-YP1))
      NSTEP  = 1 + INT(1.5*RAYLEN/DSTEP + 0.5)                            30.70
!
      DO 10 JJ = 0, NSTEP                                                 30.70
        X3  = XP1 + REAL(JJ)*(XP2-XP1)/REAL(NSTEP)
        Y3  = YP1 + REAL(JJ)*(YP2-YP1)/REAL(NSTEP)
        IF (BOTDEP) THEN
          D3   = SVALQI (X3, Y3, 1, BOTLEV, 1, 0, 0)                      30.70
        ELSE
          D3   = SVALQI (X3, Y3, 1, BOTLEV, 1, 0, 0) + WLEV               30.70
          IF (LEDS(7).GE.2)                                               30.70
     &    D3 = D3 + SVALQI (X3, Y3, 7, WATLEV, 1, 0, 0)                   30.70
        ENDIF
        IF (ITEST.GE.160) WRITE (PRTEST, 14) X3+XOFFS, Y3+YOFFS, D3       30.70
  14    FORMAT (' SIRAY, scan point', 2(1X,F8.0), 1X, F8.2)               30.70
        IF (ABS(D3-DP).LT.DIFDEP) THEN                                    10.20
          DIFDEP=ABS(D3-DP)                                               10.20
          JDMINMAX=JJ                                                     10.20
        ENDIF
        IF (JJ.GT.0) THEN
          IF ((DP-D2)*(DP-D3).LE.0) THEN                                  40.03
            IF (EQREAL(D2,D3)) THEN                                       30.72
              XX = X2
              YY = Y2
            ELSE
              XX = X2+(X3-X2)*(D2-DP)/(D2-D3)
              YY = Y2+(Y3-Y2)*(D2-DP)/(D2-D3)
            ENDIF
            SIRAY = 1
            GOTO 20
          ENDIF                                                           40.03
        ENDIF                                                             40.03
        X2 = X3
        Y2 = Y3
        D2 = D3
   10 CONTINUE
!
!     exact depth not found, take closest value:
!
      X3 = XP1 + REAL(JDMINMAX)*(XP2-XP1)/REAL(NSTEP)                     10.20
      Y3 = YP1 + REAL(JDMINMAX)*(YP2-YP1)/REAL(NSTEP)                     10.20
      XX = X3                                                             10.20
      YY = Y3                                                             10.20
!
  20  IF (ITEST.GE.140) WRITE (PRTEST, 24) XX+XOFFS, YY+YOFFS             30.70
  24  FORMAT (' SIRAY, result ', 2(1X,F8.0))                              30.70
      RETURN
! * end of function SIRAY *
      END
!************************************************************************
!                                                                       *
      SUBROUTINE SWNMPS (PSNAME, PSTYPE, MIP, IERR)                       40.31
!                                                                       *
!************************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE OUTP_DATA                                                       40.31
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
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!       Oct. 1996, ver. 30.50: new subr.
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!       Read name of set of output points; get type and number of
!       points in the set
!
!  3. METHOD
!
!
!  4. PARAMETERLIST
!
!       PSNAME   char   output   name
!       PSTYPE   char   output   type
!       MIP      int    output   number of points
!
!  5. SUBROUTINES CALLING
!
!       SPREOQ
!
!  6. SUBROUTINES USED
!
!       INCSTR (Ocean Pack)
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
!
! 10. SOURCE TEXT
!
      INTEGER   MIP                                                       40.31
      CHARACTER PSNAME *(*), PSTYPE *1                                    40.31
      TYPE(OPSDAT), POINTER :: CUOPS                                       40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SWNMPS')
!
      IERR = 0
      CALL INCSTR ('SNAME', PSNAME, 'STA', 'BOTTGRID')
      IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
      CUOPS => FOPS                                                       40.31
      DO                                                                  40.31
        IF (CUOPS%PSNAME.EQ.PSNAME) EXIT                                  40.31
        IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) THEN                          40.31
           CALL MSGERR(2, 'Set of output locations is not known')         40.31
           GOTO 900                                                       40.31
        END IF                                                            40.31
        CUOPS => CUOPS%NEXTOPS                                            40.31
      END DO                                                              40.31
      PSTYPE = CUOPS%PSTYPE                                               40.31
      IF (PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H') THEN                          40.31
         MIP = CUOPS%OPI(1) * CUOPS%OPI(2)                                40.31
!        get direction of frame in case of coordinates plotting
         ALPQ = CUOPS%OPR(5)                                              40.31
      ELSE
         MIP  = CUOPS%MIP                                                 40.31
         ALPQ = 0.                                                        40.31
      ENDIF
 800  IF (ITEST.GE.100) WRITE (PRTEST, 802) PSNAME, PSTYPE, MIP
 802  FORMAT (' exit SWNMPS, name:', A8, ' type:', A1,
     &        '  num of p:', I5)
      RETURN
 900  PSTYPE = ' '
      MIP    = 0
      IERR   = 1
      RETURN
      END
!************************************************************************
!                                                                       *
      SUBROUTINE SVARTP (IVTYPE)
!                                                                       *
!************************************************************************
!
      USE SWCOMM1                                                         40.41
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
!     32.02: Roeland Ris & Cor van der Schelde
!     40.03: Nico Booij
!     40.04: Annette Kieftenburg
!
!  1. Updates
!
!     10.09, Aug. 94: output quantity RPER added
!     20.61, Sep. 95: quantities TM02 and FWID added
!     20.67, Dec. 95: FWID renamed FSPR (freq. spread)
!     32.02, Feb. 98: 1D-version introduced
!     40.00, Apr. 98: subr simplified using new array OVKEYW
!     40.03, Sep. 00: inconsistency with manual corrected
!
!  2. Purpose
!
!     Converting keyword into integer
!
!  3. Method
!
!     This subroutine determines an integer value indicating the
!     required output variable from the keyword denoting the same
!     for storage in array with output requests.
!
!  4. PARAMETERLIST
!
!     IVTYPE  INT    output   type number output variable
!
!  5. Subroutines calling
!
!     SPROUT
!
!  6. Subroutines used
!
!     KEYWIS (Ocean Pack)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     -----------------------------------------------------------------
!     If the keyword is equal to given string, then
!         IVTYPE is given integer value
!     -----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL KEYWIS
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SVARTP')
!
      IVTYPE  =  0
!
      CALL INKEYW ('STA', 'ZZZZ')
!     check if given keyword corresponds to output quantity
      DO IVT = NMOVAR, 1, -1
!       loop in reverse order to check more specific names first
!       e.g. HSWE before HS
        IF (KEYWIS (OVKEYW(IVT))) THEN                                    40.00
          IVTYPE = IVT                                                    40.00
          GOTO 40
        ENDIF
      ENDDO
!     aliases:
      IF (KEYWIS ('PPER')) IVTYPE = 12                                    40.00
      IF (KEYWIS ('RPER')) IVTYPE = 28                                    40.00
      IF (KEYWIS ( 'DTM')) IVTYPE = 31                                    40.00
      IF (KEYWIS ('FWID')) IVTYPE = 33                                    40.00
!     keyword OUTPUT means that output times will be entered              40.00
      IF (KEYWIS ('OUT')) IVTYPE = 98                                     40.03
!     keyword ZZZZ means end of list of output quantities                 40.00
      IF (KEYWIS ('ZZZZ')) IVTYPE = 99
!
      IF (IVTYPE .EQ. 0) CALL WRNKEY
!
  40  RETURN
!     end of subroutine SVARTP *
      END
!************************************************************************
!                                                                       *
      SUBROUTINE SWBOUN ( XCGRID, YCGRID, KGRPNT, XYTST, KGRBND )         40.31
!                                                                       *
!************************************************************************
!
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE M_BNDSPEC                                                       40.31
      USE SwanGriddata                                                    40.80
      USE SwanGridobjects                                                 40.80
      USE SwanCompdata                                                    40.80
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
!     30.73: Nico Booij
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     34.01: Jeroen Adema
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.05: Ekaterini E. Kriezi
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!     40.92: Marcel Zijlema
!     41.14: Nico Booij
!
!  1. Updates
!
!     30.73, Nov. 97: New subroutine, replacing code in subr. SWREAD (file SWANPRE1)
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description several arrays
!     30.81, Nov. 98: Adjustment for 1-D case of new boundary conditions
!     34.01, Feb. 99: Introducing STPNOW
!     30.81, Apr. 99: Prevent negative powers for cosine directional spreading (DSPR);
!                     prevented DSPR > 360 and DSPR < 0 (except for exception value).
!     30.82, July 99: Used EQREAL for real equality comparisons
!     40.05, Aug  00: WW3 boundary nesting command, in Swan nesting option
!                     adding of a new option (same as WW3 command)
!     40.03, Sep. 00: inconsistency with manual corrected
!     40.02, Oct. 00: WWIII added as keyword (will appear in the manual)
!     40.13, Nov. 01: determination of side corrected (iside=3)
!     40.31, Nov. 03: removing POOL-mechanism, reconsideration of this subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Jun. 07: extension to unstructured grids
!     40.92, Jun. 08: changes with respect to boundary polygons
!     41.14, Jul. 10: call SwanBndStruc added
!
!  2. Purpose
!
!     Reading and processing BOUNDARY command
!
!  3. Method
!
!
!  4. Argument variables
!
! i   XCGRID: Coordinates of computational grid in x-direction            30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.82
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.82
!
!     KGRBND  int   inp    grid indices of boundary points                40.31
!     KGRPNT  int   inp    indirect addresses of grid points
!     XYTST   int   inp    ix, iy of test points
!
      INTEGER KGRPNT(MXC,MYC)
      INTEGER XYTST(*),  KGRBND(*)
!
!  5. Parameter variables
!
!
!  6. Local variables
!
      INTEGER   IENT,KOUNTR,IX1,IY1,IX2,IY2
      INTEGER   MM,IX,IY,ISIDM,ISIDE,KC,KC2,KC1,IX3,IY3,MP
      INTEGER   IP,II,NBSPSS,NFSEQ,IKO,IKO2,IBSPC1,IBSPC2
      INTEGER   VM                                                        40.80

      INTEGER, DIMENSION(:), ALLOCATABLE :: IARR1, IARR2                  40.92

      REAL      CRDP, CRDM, SOMX, SOMY
      REAL      XP,YP,XC,YC,RR,DIRSI,COSDIR,SINDIR,DIRSID,DIRREF
      REAL      RLEN1,RDIST,RLEN2,XC1,YC1,XC2,YC2,W1

      LOGICAL   KEYWIS, LOCGRI, CCW, BPARF, BOUNPT,DONALL
      LOGICAL   LFRST1, LFRST2, LFRST3                                    40.31

      INTEGER   NUMP

      LOGICAL, SAVE :: LBFILS = .FALSE.                                   40.31
      LOGICAL, SAVE :: LBS    = .FALSE.                                   40.31
      LOGICAL, SAVE :: LBGP   = .FALSE.                                   40.31

      TYPE(BSPCDAT), POINTER :: BFLTMP                                    40.31
      TYPE(BSPCDAT), SAVE, POINTER :: CUBFL                               40.31

      TYPE(BSDAT), POINTER :: BSTMP                                       40.31
      TYPE(BSDAT), SAVE, POINTER :: CUBS                                  40.31

      TYPE(BGPDAT), POINTER :: BGPTMP                                     40.31

      TYPE XYPT                                                           40.31
        INTEGER             :: JX, JY
        TYPE(XYPT), POINTER :: NEXTXY
      END TYPE XYPT

      TYPE(XYPT), TARGET  :: FRST                                         40.31
      TYPE(XYPT), POINTER :: CURR, TMP                                    40.31

      CHARACTER(80) :: MSGSTR                                             40.80

      TYPE(verttype), DIMENSION(:), POINTER :: vert                       40.80
      TYPE(facetype), DIMENSION(:), POINTER :: face                       40.80
!
!  8. Subroutines used
!
!       Ocean Pack command reading routines
!       BOUNPT
!
      LOGICAL STPNOW                                                      34.01
      LOGICAL EQREAL
!
!  9. Subroutines calling
!
!       SWREAD
!
!  10. Error messages
!
!
!  11. Remarks
!
!       data concerning boundary files are stored in array BFILES
!       see subr BCFILE for details
!
!  12. Structure
!
!       -----------------------------------------------------------------
!       Read keyword
!       Case keyword =
!       'SHAPE': Read spectral shape parameters
!       'WAMN': Read filename
!               Open WAM nesting file
!               Put file characteristics into array BFILES
!       'WW3N': Read filename
!               Open WW3 nesting file
!               Put file characteristics into array BFILES
!       'NEST': Read filename
!               Call BCFILE to obtain file characteristics
!       'SIDE': read side of boundary
!       'SEG':  read set of points on boundary of comp. grid
!               Read keyword
!               If keyword is 'UNIF'
!               Then Read keyword
!                    If keyword is 'PAR'
!                    Then read integral wave parameters
!                         Call SSHAPE to generate spectrum
!                         put spectrum into array BSPECS
!                    Else {keyword is 'FILE'}
!                         Read filename
!                         Call BCFILE to obtain file characteristics
!               Else {keyword is 'VAR'}
!                    Read keyword
!                    If keyword is 'PAR'
!                    Then Repeat until list is exhausted
!                             Read length and integral wave parameters
!                             Call SSHAPE to generate spectrum
!                             put spectrum into array BSPECS
!                    Else {keyword is 'FILE'}
!                         Repeat until list is exhausted
!                             Read filename
!                             Call BCFILE to obtain file characteristics
!       -----------------------------------------------------------------
!
! 13. Source text

      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SWBOUN')
!
!     point to vertex and face objects
!
      vert => gridobject%vert_grid                                        40.80
      face => gridobject%face_grid                                        40.80
!
      IF (OPTG.EQ.5) THEN                                                 40.80
!
!        in case of unstructured grid, make list of boundary points
!        in ascending order
!
         CALL SwanBpntlist                                                40.80
      ELSE                                                                41.14
!
!        generate output curves BOUNDARY and BOUND_** for structured grids
!
         CALL SwanBndStruc ( XCGRID, YCGRID )                             41.14
      ENDIF                                                               40.80
!
      CALL INKEYW ('REQ',' ')
      IF (KEYWIS ('SHAP')) THEN
!
!           specification of the spectral shape
!
! =========================================================================
!
!                      |  JONswap  [gamma]  |
!                      |                    |    | -> PEAK |
!  BOUNdspec  SHAPe   <   PM                 >  <           >   &
!                      |                    |    | MEAN    |
!                      |  GAUSs  [sigfr]    |
!                      |                    |
!                      |  BIN               |
!
!                     | DEGRees   |
!             DSPR   <             >
!                     | -> POWer  |
!
! =========================================================================
!
        CALL INKEYW ('STA', 'JON')
        IF (KEYWIS ('JON')) THEN
          FSHAPE = 2
          CALL INREAL ('GAMMA', PSHAPE(1), 'STA', 3.3)                    40.00
        ELSE IF (KEYWIS ('BIN')) THEN
          FSHAPE = 3
        ELSE IF (KEYWIS ('PM')) THEN
          FSHAPE = 1
        ELSE IF (KEYWIS ('GAUS')) THEN
          FSHAPE = 4
          CALL INREAL ('SIGFR', SIGMAG, 'STA', 0.01)
!         convert from Hz to rad/s:
          PSHAPE(2) = PI2 * SIGMAG                                        40.00
        ENDIF
!       PEAK or MEAN frequency
        CALL INKEYW ('STA', ' ')
        IF (KEYWIS('MEAN')) THEN
          FSHAPE = -FSHAPE
        ELSE
          CALL IGNORE ('PEAK')
        ENDIF
!       directional distribution given by DEGR or by POWER
        CALL IGNORE ('DSPR')
        CALL INKEYW ('STA', 'POW')
        IF (KEYWIS('DEGR')) THEN
          DSHAPE = 1
        ELSE
          CALL IGNORE ('POW')
          DSHAPE = 2
        ENDIF
        IF (ITEST.GE.30) WRITE (PRINTF,6100) FSHAPE, DSHAPE
 6100   FORMAT (' Shape of inc. spectrum, Freq:', I2, ' ; Dir:', I2)
!
      ELSE IF (KEYWIS ('WAMN')) THEN
!
!       SWAN in WAM nesting                                               30.04
!
!      =======================================================================
!
!                                                   |-> CRAY |
!                                    | UNFormatted <          > |
!                                    |              | WKstat |  |
!                                    |                          |
!       BOUNdnest2  WAMNest 'fname' <                            > [xgc] [ygc] [lwdate]
!                                    |                          |
!                                    | FREE                     |
!
!      =======================================================================
!
        IF (MXC .LE. 0 .AND. OPTG.NE.5) THEN                              40.80
          CALL MSGERR(3, ' command CGRID must precede this command')      40.80
          GOTO 900
        ENDIF
        IF (MCGRD .LE. 1 .AND. nverts .LE. 0) THEN                        40.80
          CALL MSGERR(3,
     &    ' command READ BOT or READ UNSTRUC must precede this command')  40.80
          GOTO 900
        ENDIF
!
        IF (OPTG.EQ.5) THEN                                               40.80
           CALL MSGERR(2,
     &               ' WAM b.c. are not supported in unstructured grid')  40.80
           GOTO 900                                                       40.80
        ENDIF                                                             40.80
!
        NBFILS = NBFILS + 1
        ALLOCATE(BFLTMP)                                                  40.31
        CALL INCSTR ('FNAME',FILENM,'REQ', ' ')
        CALL BCWAMN (FILENM, 'NEST', BFLTMP, LBGP,                        40.31
     &               XCGRID, YCGRID, KGRPNT, XYTST)                       40.31
        IF (STPNOW()) RETURN                                              34.01
        NULLIFY(BFLTMP%NEXTBSPC)                                          40.31
        IF ( .NOT.LBFILS ) THEN                                           40.31
           FBNDFIL = BFLTMP                                               40.31
           CUBFL => FBNDFIL                                               40.31
           LBFILS = .TRUE.                                                40.31
        ELSE                                                              40.31
           CUBFL%NEXTBSPC => BFLTMP                                       40.31
           CUBFL => BFLTMP                                                40.31
        END IF                                                            40.31
!
      ELSE IF (KEYWIS('WW3').OR.KEYWIS('WWIII')) THEN                     40.02
!
!       SWAN in WaveWatch nesting                                         40.05
!
!      =======================================================================
!                                 | -> CLOS |
!       BOUNdnest2  WWIII 'fname'  <          > [xgc] [ygc]               40.02
!                                 |  OPEN   |
!      =======================================================================
!
        IF (MXC .LE. 0 .AND. OPTG.NE.5) THEN                              40.80
          CALL MSGERR(3, ' command CGRID must precede this command')      40.80
          GOTO 900
        ENDIF
        IF (MCGRD .LE. 1 .AND. nverts .LE. 0) THEN                        40.80
          CALL MSGERR(3,
     &    ' command READ BOT or READ UNSTRUC must precede this command')  40.80
          GOTO 900
        ENDIF
!
        IF (OPTG.EQ.5) THEN                                               40.80
           CALL MSGERR(2,
     &             ' WWIII b.c. are not supported in unstructured grid')  40.80
           GOTO 900                                                       40.80
        ENDIF                                                             40.80
!
        NBFILS = NBFILS + 1
        ALLOCATE(BFLTMP)                                                  40.31
        CALL INCSTR ('FNAME',FILENM,'STA', 'nest.ww3')

!       if keyword is OPEN (boundary is not a closed contour) then
!          DONALL is TRUE  and the nesting boundary remain open
!       else (default case)
!          DONALL is FALSE  and boundary is close and interpolation
!          between the last and the first point will be done

        CALL INKEYW ('STA', 'CLOS')                                       40.05
        IF (KEYWIS('OPEN')) THEN                                          40.05
           DONALL = .TRUE.                                                40.05
        ELSE IF (KEYWIS('CLOS')) THEN                                     40.05
           DONALL = .FALSE.                                               40.05
        ELSE                                                              40.05
           CALL WRNKEY                                                    40.05
        ENDIF                                                             40.05
        CALL BCWW3N (FILENM, 'NEST', BFLTMP, LBGP,                        40.31 40.05
     &               XCGRID, YCGRID, KGRPNT, XYTST, KGRBND,               40.31 40.05
     &               DONALL)                                              40.31 40.05
        IF (STPNOW()) RETURN
        NULLIFY(BFLTMP%NEXTBSPC)                                          40.31
        IF ( .NOT.LBFILS ) THEN                                           40.31
           FBNDFIL = BFLTMP                                               40.31
           CUBFL => FBNDFIL                                               40.31
           LBFILS = .TRUE.                                                40.31
        ELSE                                                              40.31
           CUBFL%NEXTBSPC => BFLTMP                                       40.31
           CUBFL => BFLTMP                                                40.31
        END IF                                                            40.31
!
      ELSE IF (KEYWIS ('NE')) THEN
!
!       Nesting SWAN model in larger SWAN model
! ==========================================
!                                | -> CLOS |                              40.05
!     BOUNdnest1  NEST 'fname'  <           >
!                                |  OPEN   |                              40.05
! ==========================================
!
        IF (MXC .LE. 0 .AND. OPTG.NE.5) THEN                              40.80
          CALL MSGERR(3, ' command CGRID must precede this command')      40.80
          GOTO 900
        ENDIF
        IF (MCGRD .LE. 1 .AND. nverts .LE. 0) THEN                        40.80
          CALL MSGERR(3,
     &    ' command READ BOT or READ UNSTRUC must precede this command')  40.80
          GOTO 900
        ENDIF
!
        NBFILS = NBFILS + 1
        ALLOCATE(BFLTMP)                                                  40.31
        CALL INCSTR ('FNAME',FILENM,'REQ', ' ')

!       if keyword is OPEN then
!          DONALL is TRUE  and the nesting boundary remain open
!       else (default case)
!          DONALL is FALSE  and boundary is close and interpolation between
!          the last and the first point will be done

        CALL INKEYW ('STA', 'CLOS')                                       40.05
        IF (KEYWIS('OPEN')) THEN                                          40.05
           DONALL = .TRUE.                                                40.05
        ELSE IF (KEYWIS('CLOS')) THEN                                     40.05
           DONALL = .FALSE.                                               40.05
        ELSE                                                              40.05
           CALL WRNKEY                                                    40.05
        ENDIF                                                             40.05

        CALL BCFILE (FILENM, 'NEST', BFLTMP, LBGP,                        40.31
     &               XCGRID, YCGRID, KGRPNT, XYTST,  KGRBND,              40.31
     &               DONALL)                                              40.05
        IF (STPNOW()) RETURN                                              34.01
        NULLIFY(BFLTMP%NEXTBSPC)                                          40.31
        IF ( .NOT.LBFILS ) THEN                                           40.31
           FBNDFIL = BFLTMP                                               40.31
           CUBFL => FBNDFIL                                               40.31
           LBFILS = .TRUE.                                                40.31
        ELSE                                                              40.31
           CUBFL%NEXTBSPC => BFLTMP                                       40.31
           CUBFL => BFLTMP                                                40.31
        END IF                                                            40.31
!
      ELSE
!
!       parametric or file boundary condition
!
!      ======================================================================
!
!                             | North |
!                             | NW    |
!                             | West  |
!                             | SW    |          | -> CCW     |
!                 | -> SIDE  <  South  > | [k]  <              >   |
!                 |           | SE    |          | CLOCKWise  |    |
!                 |           | East  |                            |
!                 |           | NE    |                            |
!       BOUNdary <                                                  >    &
!                 |           | -> XY  < [x] [y] >           |     |
!                 | SEGment  <                                >    |
!                             |    IJ  < [i] [j] > | < [k] > |
!
!
!                          |  PAR  [hs] [per] [dir] [dd]  |
!            |  UNIForm   <                                >             |
!            |             |  FILE  'fname'  [seq]        |              |
!           <                                                             >
!            |             |  PAR  < [len] [hs] [per] [dir] [dd] >  |    |
!            |  VARiable  <                                          >   |
!                          |  FILE < [len] 'fname' [seq] >          |
!
!      ======================================================================
!
        IF (MXC .LE. 0 .AND. OPTG.NE.5) THEN                              40.80
          CALL MSGERR(3, ' command CGRID must precede this command')      40.80
          GOTO 900
        ENDIF
        IF (MCGRD .LE. 1 .AND. nverts .LE. 0) THEN                        40.80
          CALL MSGERR(3,
     &    ' command READ BOT or READ UNSTRUC must precede this command')  40.80
          GOTO 900
        ENDIF
!
!       first define side or segment
!
!       *** definition of boundary segment ***
!
        CALL INKEYW ('REQ',' ')
        IF (KEYWIS ('STAT')) THEN
          CALL MSGERR (1, 'keyword STAT ignored')
          CALL INKEYW ('REQ',' ')
        ENDIF
        KOUNTR  = 0
        FRST%JX = 0                                                       40.31
        FRST%JY = 0                                                       40.31
        NULLIFY(FRST%NEXTXY)                                              40.31
        CURR => FRST                                                      40.31
        IF (KEYWIS ('SEG')) THEN
          IERR = 0
          CALL INKEYW ('STA','XY')
          IF (KEYWIS('XY') .OR. KEYWIS ('LOC')) THEN
            LOCGRI = .TRUE.
          ELSE IF (KEYWIS('IJ') .OR. KEYWIS ('GRI')) THEN
            LOCGRI = .FALSE.
          ELSE
            CALL WRNKEY
          ENDIF
          IX1 = 1
          IY1 = 1
          LFRST1 = .TRUE.                                                 40.31
!         loop over points describing the segment
          DO
            IF (LOCGRI) THEN
              CALL READXY ('XP','YP',XP,YP, 'REP', -1.E10, -1.E10)
              IF (XP.LT.-.9E10) GOTO 42
              IF (OPTG.NE.5) THEN                                         40.80
!             --- structured grid
!
                 CALL CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID, YCGRID,
     &                        KGRBND)                                     40.00
                 IX2 = NINT(XC) + 1
                 IY2 = NINT(YC) + 1                                       40.00
                 IF (.NOT.BOUNPT(IX2,IY2,KGRPNT)) THEN
                   CALL MSGERR (2, 'invalid boundary point')
                   WRITE (PRTEST, 38) XP+XOFFS, YP+YOFFS, XC, YC,         40.00
     &                                IX2, IY2                            40.00
  38               FORMAT (' segment point ', 2F10.2,
     &                     ' grid ', 2F8.2, 2I4)
                 ENDIF
              ELSE                                                        40.80
!             --- unstructured grid
!
                 CALL SwanFindPoint ( XP, YP, IX2 )                       40.80
                 IF ( IX2.LT.0 ) THEN                                     40.80
                    WRITE (MSGSTR, '(A,F12.4,A,F12.4,A)')                 40.80
     &                        ' Boundary point (',XP+XOFFS,',',YP+YOFFS,  40.80
     &                                ') not part of computational grid'  40.80
                    CALL MSGERR( 2, TRIM(MSGSTR) )                        40.80
                 ENDIF                                                    40.80
                 IF ( vert(IX2)%atti(VMARKER) /= 1 ) THEN                 40.80
                    WRITE (MSGSTR, '(A,F12.4,A,F12.4,A)')                 40.80
     &                                ' Vertex (',XP+XOFFS,',',YP+YOFFS,  40.80
     &                                 ') is not a valid boundary point'  40.80
                    CALL MSGERR( 2, TRIM(MSGSTR) )                        40.80
                 ENDIF                                                    40.80
              ENDIF                                                       40.80
            ELSE
              CALL ININTG ('I' , IX2, 'REP', -1)                          40.03
              IF (IX2 .LT. 0) GOTO 42                                     40.00
              IF (OPTG.NE.5) THEN                                         40.80
                 CALL ININTG ('J' , IY2, 'REQ',  0)                       40.03
                 IX2 = IX2 + 1                                            40.00
                 IY2 = IY2 + 1
              ELSE                                                        40.80
                 IF (IX2.LE.0 .OR. IX2.GT.nverts) THEN                    40.80
                    WRITE (MSGSTR,'(I4,A)') IX2,                          40.80
     &                                    ' is not a valid vertex index'  40.80
                    CALL MSGERR( 2, TRIM(MSGSTR) )                        40.80
                 ELSEIF ( vert(IX2)%atti(VMARKER) /= 1 ) THEN             40.80
                    WRITE (MSGSTR,'(A,I4,A)') ' Vertex with index ',IX2,  40.80
     &                                  ' is not a valid boundary point'  40.80
                    CALL MSGERR( 2, TRIM(MSGSTR) )                        40.80
                 ENDIF                                                    40.80
              ENDIF                                                       40.80
            ENDIF
            IF (ITEST.GE.80 .AND. OPTG.NE.5) WRITE (PRTEST, 38)           40.00
     &                XCGRID(IX2,IY2)+XOFFS,                              40.00
     &                YCGRID(IX2,IY2)+YOFFS, XC, YC, IX2-1, IY2-1         40.00
!
!           --- generate intermediate points on the segment
            IF ( OPTG.NE.5 ) THEN                                         40.80
!           --- structured grid                                           40.80
!
               IF (IX2 .GT. 0 .AND. IX2 .LE. MXC .AND.
     &             IY2 .GT. 0 .AND. IY2 .LE. MYC) THEN
                 IF (LFRST1) THEN
                   MM = 1
                   LFRST1 = .FALSE.
                 ELSE
                   MM = MAX (ABS(IX2-IX1), ABS(IY2-IY1))
                 ENDIF
                 DO IP = 1, MM
                   RR = REAL(IP) / REAL(MM)
!
                   IF (.NOT. ONED) THEN                                   30.81
                     IX = IX1 + NINT(RR*REAL(IX2-IX1))
                     IY = IY1 + NINT(RR*REAL(IY2-IY1))
                   ELSE                                                   30.81
                     IX = IX1 + NINT(RR*REAL(IX2-IX1))                    30.81
                     IY = IY1                                             30.81
                   END IF                                                 30.81
!
                   IF (ITEST.GE.80) WRITE (PRTEST, *) ' b. point ',
     &                   RR, IX, IY
                   IF (KGRPNT(IX,IY) .GT. 1) THEN
                     KOUNTR = KOUNTR + 1
                     ALLOCATE(TMP)                                        40.31
                     TMP%JX = IX                                          40.31
                     TMP%JY = IY                                          40.31
                     NULLIFY(TMP%NEXTXY)                                  40.31
                     CURR%NEXTXY => TMP                                   40.31
                     CURR => TMP                                          40.31
                   ENDIF
                 ENDDO
               ELSE
                 MSGSTR =''                                               41.14
                 write (MSGSTR, 117) IX2-1, IY2-1                         41.14
 117             format ('(',2I5, ') is outside computational grid')      41.14
                 CALL MSGERR (2, MSGSTR)
               ENDIF
               IY1 = IY2
            ELSE                                                          40.80
!           --- unstructured grid                                         40.80
!
               IF (LFRST1) THEN                                           40.80
                  IXB1 = vert(IX2)%atti(BINDX)                            40.80
                  JBG  = vert(IX2)%atti(BPOL)                             40.92
                  DET  = 1.                                               40.92
                  LFRST1 = .FALSE.                                        40.80
               ELSE                                                       40.80
                  IXB1 = vert(IX1)%atti(BINDX)                            40.80
                  JBG  = vert(IX1)%atti(BPOL)                             40.92
!
!                 1) the wave spectrum along the given segment can be
!                    imposed in counterclockwise or clockwise direction
!                 2) content of array blist is ordered in counterclockwise
!                    manner for sea/mainland boundary (JBG=1) and
!                    clockwise for island boundary (JBG>1)
!                 3) therefore, determine orientation by means of the
!                    determinant of two endpoints of the given segment
!                    and an arbitrary point inside domain
!
!                 first endpoint of segment
                  X1 = vert(IX1)%attr(VERTX)                              40.80
                  Y1 = vert(IX1)%attr(VERTY)                              40.80
!
!                 second endpoint of segment
                  X2 = vert(IX2)%attr(VERTX)                              40.80
                  Y2 = vert(IX2)%attr(VERTY)                              40.80
!
!                 an arbitrary internal point                             40.92
                  DX=0.1*mingsiz
                  DY=0.1*mingsiz
                  DO IP=1,4
                     X3 = X1 - DX
                     Y3 = Y1 - DY
                     CALL SwanFindPoint ( X3, Y3, IX )
                     IF ( JBG>1 .AND. IX.LT.0 ) THEN
                        X3 = X1 + DX
                        Y3 = Y1 + DY
                        EXIT
                     ELSEIF ( JBG==1 .AND. IX.GT.0 ) THEN
                        EXIT
                     ENDIF
                     IF ( MOD(IP,2).EQ.0 ) THEN
                        DX = -DX
                     ELSE
                        DY = -DY
                     ENDIF
                  ENDDO
!
                  DET= (Y3-Y1)*(X2-X1)-(Y2-Y1)*(X3-X1)                    40.80
                  IF (DET.GT.0.) THEN                                     40.80
!                    take next boundary point in counterclockwise
!                    direction
                     IXB1 = MOD(IXB1,nbpt(JBG))+1                         40.92 40.80
                  ELSE                                                    40.80
!                    take next boundary point in clockwise direction
                     IXB1 = nbpt(JBG)-MOD(nbpt(JBG)+1-IXB1,nbpt(JBG))     40.92 40.80
                  ENDIF                                                   40.80
               ENDIF                                                      40.80
               IXB2 = vert(IX2)%atti(BINDX)                               40.80
!
!              determine order of counting
               IF (IXB1.GT.IXB2 ) THEN                                    40.80
                  IF (DET.LT.0.) THEN                                     40.92
                     IXI  = -1                                            40.80
                  ELSE                                                    40.92
                     IXI  = 1                                             40.92
                     IXB2 = IXB2+nbpt(JBG)                                40.92
                  ENDIF                                                   40.92
               ELSE                                                       40.80
                  IF (DET.GT.0.) THEN                                     40.92
                     IXI  = 1                                             40.80
                  ELSE                                                    40.92
                     IXI  = -1                                            40.92
                     IXB1 = IXB1+nbpt(JBG)                                40.92
                  ENDIF                                                   40.92
               ENDIF                                                      40.80
!
               DO IPP = IXB1, IXB2, IXI                                   40.92 40.80
                  IP = MOD(IPP,nbpt(JBG))                                 40.92
                  IF (IP.EQ.0) IP = nbpt(JBG)                             40.92
                  KOUNTR = KOUNTR + 1                                     40.80
                  IX = blist(IP,JBG)                                      40.92 40.80
                  vert(IX)%atti(VBC) = 1                                  40.80
                  ALLOCATE(TMP)                                           40.80
                  TMP%JX = IX                                             40.80
                  NULLIFY(TMP%NEXTXY)                                     40.80
                  CURR%NEXTXY => TMP                                      40.80
                  CURR => TMP                                             40.80
               ENDDO                                                      40.80
            ENDIF                                                         40.80
            IX1 = IX2
          END DO
 42       IF (KOUNTR.EQ.0)
     &                CALL MSGERR(1,'No points on the boundaries found')
          IF (LFRST1) CALL MSGERR (1,
     &        'At least two points needed for a segment')                 40.81
        ELSE
!         boundary condition on one side of the computational grid
          IF (OPTG.EQ.3) THEN                                             40.80 40.31
             CALL MSGERR(2,                                               40.80 40.31
     &          ' keyword SIDE should not be used for curvilinear grid')  40.80 40.31
          END IF                                                          40.31
          CALL IGNORE ('SIDE')
!         *** specification of side for which boundary   ***
!         *** condition is given                         ***
          IF (OPTG.NE.5) THEN                                             40.80
             CALL INKEYW ('REQ',' ')
             IF (KEYWIS ('NW')) THEN
               DIRSI = 45.
             ELSE IF (KEYWIS ('SW')) THEN
               DIRSI = 135.
             ELSE IF (KEYWIS ('SE')) THEN
               DIRSI = -135.
             ELSE IF (KEYWIS ('NE')) THEN
               DIRSI = -45.
             ELSE IF (KEYWIS ('N')) THEN
               DIRSI = 0.
             ELSE IF (KEYWIS ('W')) THEN
               DIRSI = 90.
             ELSE IF (KEYWIS ('S')) THEN
               DIRSI = 180.
             ELSE IF (KEYWIS ('E')) THEN
               DIRSI = -90.
             ELSE
               CALL WRNKEY
             ENDIF
          ELSE                                                            40.80
             CALL ININTG ('K', VM, 'REQ', 0)                              40.80
          ENDIF                                                           40.80
!
!         --- go along boundary clockwise or counterclockwise (default)
!
          CALL INKEYW ('STA', 'CCW')
          IF (KEYWIS('CLOCKW')) THEN
            CCW = .FALSE.
          ELSE
            CALL IGNORE ('CCW')
            CCW = .TRUE.
          ENDIF
!
!         select side in the chosen direction
!
          IF ( OPTG.NE.5 ) THEN                                           40.80
             CRDM   = -1.E10
             ISIDM  = 0
             IF (ONED) THEN                                               40.00
               COSDIR = COS(PI*(DNORTH+DIRSI)/180.)
               SINDIR = SIN(PI*(DNORTH+DIRSI)/180.)
               DO ISIDE = 1, 4
                 SOMX = 0.
                 SOMY = 0.
                 NUMP = 0
                 IF (ISIDE.EQ.2) THEN
                   KC = KGRPNT(MXC,1)
                   IF (KC.GT.1) THEN
                     SOMX = XCGRID(MXC,1)
                     SOMY = YCGRID(MXC,1)
                     NUMP = 1
                   ENDIF
                 ELSE IF (ISIDE.EQ.4) THEN
                   KC = KGRPNT(1,1)
                   IF (KC.GT.1) THEN
                     SOMX = XCGRID(1,1)
                     SOMY = YCGRID(1,1)
                     NUMP = 1
                   ENDIF
                 ENDIF
                 IF (NUMP.GT.0) THEN
                   CRDP = COSDIR*SOMX + SINDIR*SOMY
!                  side with largest CRDP is the one selected
                   IF (CRDP.GT.CRDM) THEN
                     CRDM = CRDP
                     ISIDM = ISIDE
                   ENDIF
                 ENDIF
               ENDDO
             ELSE                                                         40.00
               DO ISIDE = 1, 4
                 SOMX = 0.
                 SOMY = 0.
                 NUMP = 0
                 IF (ISIDE.EQ.1) THEN
                   DO IX = 1, MXC
                     KC2 = KGRPNT(IX,1)
                     IF (IX.GT.1) THEN
                       IF (KC1.GT.1 .AND. KC2.GT.1) THEN                  40.00
!                        if both grid points at ends of a step are valid, then
!                        take DX and DY into account when determining direction
                         SOMX = SOMX + XCGRID(IX,1)-XCGRID(IX-1,1)
                         SOMY = SOMY + YCGRID(IX,1)-YCGRID(IX-1,1)
                         NUMP = NUMP + 1
                       ENDIF
                     ENDIF
                     KC1 = KC2                                            40.03
                   ENDDO
                 ELSE IF (ISIDE.EQ.2) THEN
                   DO IY = 1, MYC
                     KC2 = KGRPNT(MXC,IY)
                     IF (IY.GT.1) THEN
                       IF (KC1.GT.1 .AND. KC2.GT.1) THEN                  40.00
                         SOMX = SOMX + XCGRID(MXC,IY)-XCGRID(MXC,IY-1)
                         SOMY = SOMY + YCGRID(MXC,IY)-YCGRID(MXC,IY-1)
                         NUMP = NUMP + 1
                       ENDIF
                     ENDIF
                     KC1 = KC2                                            40.03
                   ENDDO
                 ELSE IF (ISIDE.EQ.3) THEN                                40.00
                   DO IX = 1, MXC
                     KC2 = KGRPNT(IX,MYC)
                     IF (IX.GT.1) THEN
                       IF (KC1.GT.1 .AND. KC2.GT.1) THEN                  40.00
                         SOMX = SOMX + XCGRID(IX-1,MYC)-XCGRID(IX,MYC)    40.13
                         SOMY = SOMY + YCGRID(IX-1,MYC)-YCGRID(IX,MYC)    40.13
                         NUMP = NUMP + 1
                       ENDIF
                     ENDIF
                     KC1 = KC2                                            40.03
                   ENDDO
                 ELSE IF (ISIDE.EQ.4) THEN
                   DO IY = 1, MYC
                     KC2 = KGRPNT(1,IY)
                     IF (IY.GT.1) THEN
                       IF (KC1.GT.1 .AND. KC2.GT.1) THEN                  40.00
                         SOMX = SOMX + XCGRID(1,IY-1)-XCGRID(1,IY)
                         SOMY = SOMY + YCGRID(1,IY-1)-YCGRID(1,IY)
                         NUMP = NUMP + 1
                       ENDIF
                     ENDIF
                     KC1 = KC2                                            40.03
                   ENDDO
                 ENDIF
                 IF (NUMP.GT.0) THEN
                   DIRSID = ATAN2(SOMY,SOMX)
                   DIRREF = PI*(DNORTH+DIRSI)/180.
                   IF (CVLEFT) THEN
                     CRDP = COS(DIRSID - 0.5*PI - DIRREF)
                   ELSE
                     CRDP = COS(DIRSID + 0.5*PI - DIRREF)
                   ENDIF
!                  side with largest CRDP is the one selected
                   IF (CRDP.GT.CRDM) THEN
                     CRDM = CRDP
                     ISIDM = ISIDE
                   ENDIF
                 ENDIF
                 IF (ITEST.GE.60) WRITE (PRTEST, 151) ISIDE, NUMP,
     &           SOMX, SOMY, DIRSID*180/PI, DIRREF*180/PI, CRDP, CVLEFT   40.13
 151             FORMAT (' side ', 2I4, 2(1X,E11.4), 2(1X,F5.0), 2X,
     &                   F6.3, 2X, L1)                                    40.13
               ENDDO
             ENDIF                                                        40.00
             IF (ISIDM.EQ.0) THEN
               CALL MSGERR (2, 'No open boundary found')
             ENDIF
!
  90         IF (ISIDM.EQ.1) THEN
               IX1 = 1
               IY1 = 1
               IX2 = MXC
               IY2 = 1
             ELSE IF (ISIDM.EQ.2) THEN
               IX1 = MXC
               IY1 = 1
               IX2 = MXC
               IY2 = MYC
             ELSE IF (ISIDM.EQ.3) THEN
               IX1 = MXC
               IY1 = MYC
               IX2 = 1
               IY2 = MYC
             ELSE IF (ISIDM.EQ.4) THEN
               IX1 = 1
               IY1 = MYC
               IX2 = 1
               IY2 = 1
             ENDIF
             IF (.NOT.CCW .EQV. CVLEFT) THEN
!              swap end points
               IX3 = IX1
               IY3 = IY1
               IX1 = IX2
               IY1 = IY2
               IX2 = IX3
               IY2 = IY3
             ENDIF
             IF (ITEST.GE.50) WRITE (PRINTF, 112) ISIDM,
     &       IX1-1, IY1-1, XCGRID(IX1,IY1)+XOFFS, YCGRID(IX1,IY1)+YOFFS,  40.00
     &       IX2-1, IY2-1, XCGRID(IX2,IY2)+XOFFS, YCGRID(IX2,IY2)+YOFFS   40.00
 112         FORMAT (' Selected side:', I2, ' from ', 2I4, 2F9.0,         40.00
     &       ' to ', 2I4, 2F9.0)                                          40.00
             MP = MAX(ABS(IX2-IX1),ABS(IY2-IY1))
             DO IP = 0, MP
               IF (MP.EQ.0) THEN                                          40.00
                 RR = 0.
               ELSE
                 RR = REAL(IP) / REAL(MP)
               ENDIF
               IX = IX1 + NINT(RR*REAL(IX2-IX1))
               IY = IY1 + NINT(RR*REAL(IY2-IY1))
               IF (KGRPNT(IX,IY) .GT. 1) THEN
                 KOUNTR = KOUNTR + 1
                 ALLOCATE(TMP)                                            40.31
                 TMP%JX = IX                                              40.31
                 TMP%JY = IY                                              40.31
                 NULLIFY(TMP%NEXTXY)                                      40.31
                 CURR%NEXTXY => TMP                                       40.31
                 CURR => TMP                                              40.31
               ENDIF
             ENDDO
          ELSE                                                            40.80
            ! unstructured grid
            !
            DO JBG = 1, nbpol                                             40.92
               !
               ! first boundary polyogon is assumed an outer one
               ! (sea/mainland boundary) and hence, content of blist
               ! is ordered in counterclockwise manner
               !
               IF ( JBG==1 .EQV. CCW ) THEN                               40.92 40.80
                  IXB1 = 1                                                40.80
                  IXB2 = nbpt(JBG)                                        40.92 40.80
                  IXI  = 1                                                40.80
               ELSE                                                       40.80
                  IXB1 = nbpt(JBG)                                        40.92 40.80
                  IXB2 = 1                                                40.80
                  IXI  = -1                                               40.80
               ENDIF                                                      40.80
               !
               ALLOCATE(IARR1(SUM(nbpt)))
               K = 0
               DO IP = IXB1, IXB2, IXI                                    40.92
                  IX = blist(IP,JBG)                                      40.92
                  IF ( vmark(IX) == VM ) THEN                             40.92
                     K = K+1                                              40.92
                     IARR1(K) = IP                                        40.92
                  ENDIF
               ENDDO
               !
               IF ( K/=0 ) THEN                                           40.92
                  !
                  ALLOCATE(IARR2(K))                                      40.92
                  IARR2(1:K) = IARR1(1:K)                                 40.92
                  ISH = 0                                                 40.92
                  DO IPP = 2, K                                           40.92
                     IF ( IARR2(IPP)/=IARR2(IPP-1)+IXI ) THEN             40.92
                        ISH = IPP-1                                       40.92
                        EXIT                                              40.92
                     ENDIF                                                40.92
                  ENDDO                                                   40.92
                  IARR2 = CSHIFT(IARR2,ISH)                               40.92
                  !
                  DO IPP = 1, K                                           40.92
                     IP = IARR2(IPP)                                      40.92
                     IX = blist(IP,JBG)                                   40.92 40.80
                     KOUNTR = KOUNTR + 1                                  40.80
                     vert(IX)%atti(VBC) = 1                               40.80
                     ALLOCATE(TMP)                                        40.80
                     TMP%JX = IX                                          40.80
                     NULLIFY(TMP%NEXTXY)                                  40.80
                     CURR%NEXTXY => TMP                                   40.80
                     CURR => TMP                                          40.80
                  ENDDO                                                   40.80
                  DEALLOCATE(IARR2)                                       40.92
                  !
               ENDIF                                                      40.92
               DEALLOCATE(IARR1)                                          40.92
               !
            ENDDO                                                         40.92
            !
          ENDIF                                                           40.80
        ENDIF
!
!       *** boundary condition from file, 1-d or 2-d spectrum
!
        CURR => FRST%NEXTXY                                               40.31
        CALL INKEYW ('REQ',' ')
        IF (KEYWIS('UNIF') .OR. KEYWIS('CON') .OR. KEYWIS('PAR')) THEN
          CALL INKEYW('STA', 'PAR')
          IF (KEYWIS('PAR')) THEN
            CALL INREAL ('HS',  SPPARM(1), 'REQ', 0.)
            CALL INKEYW ('STA', ' ')
            CALL INREAL ('PER', SPPARM(2), 'REQ', 0.)
            CALL INREAL ('DIR', SPPARM(3), 'REQ', 0.)
            IF (DSHAPE.EQ.1) THEN
              CALL INREAL ('DD',  SPPARM(4), 'STA', 30.)
              IF ((SPPARM(4).GT.360. .OR. SPPARM(4).LT. 0.).AND.          30.81
     &            .NOT.(EQREAL(SPPARM(4),OVEXCV(16)))) THEN               30.82
                CALL MSGERR (2,'Directional spreading is less than '//    30.81
     &                        '0 or larger than 360 degrees, and no '//   30.81
     &                        'exception value')                          30.81
              END IF                                                      30.81
            ELSE
              CALL INREAL ('DD',  SPPARM(4), 'STA', 2.)
              IF (SPPARM(4).LE. 0.) THEN                                  30.81
                CALL MSGERR (2,                                           30.81
     &          'Power of cosine is less or equal to zero')               30.81
              END IF                                                      30.81
              IF (SPPARM(4)*DDIR**2/2. .GT. 1.) THEN                      40.03
                CALL MSGERR (2,                                           40.03
     &          'distribution too narrow to be represented properly')     40.03
                WRITE (PRINTF, 142) SQRT(2./SPPARM(4))*180./PI            40.03
 142            FORMAT (' Advise: choose Dtheta < ', F8.3, ' degr')       40.03
              END IF                                                      40.03
            ENDIF
            NBSPEC = NBSPEC + 1
            IF (ITEST.GE.80) WRITE (PRTEST,*) ' bound. spectr.',
     &                   NBSPEC, (SPPARM(II), II=1,4)
            NBSPSS = NBSPEC
            ALLOCATE(BSTMP)                                               40.31
            BSTMP%NBS    = NBSPEC                                         40.31
            BSTMP%FSHAPE = FSHAPE                                         40.31
            BSTMP%DSHAPE = DSHAPE                                         40.31
            BSTMP%SPPARM(1:4) = SPPARM(1:4)                               40.31
            NULLIFY(BSTMP%NEXTBS)                                         40.31
            IF ( .NOT.LBS ) THEN                                          40.31
               FBS = BSTMP                                                40.31
               CUBS => FBS                                                40.31
               LBS = .TRUE.                                               40.31
            ELSE                                                          40.31
               CUBS%NEXTBS => BSTMP                                       40.31
               CUBS => BSTMP                                              40.31
            END IF                                                        40.31
          ELSE IF (KEYWIS('FILE') .OR. KEYWIS('SPEC')) THEN
            CALL INCSTR ('FNAME',FILENM,'REQ', ' ')
!           generate new set of file data
            NBFILS = NBFILS + 1
            NBSPSS = NBSPEC
            ALLOCATE(BFLTMP)                                              40.31
            CALL BCFILE (FILENM, 'PNTS', BFLTMP, LBGP,                    40.31
     &                   XCGRID, YCGRID, KGRPNT, XYTST, KGRBND,           40.31
     &                   DONALL)                                          40.05
            IF (STPNOW()) RETURN                                          34.01
            NULLIFY(BFLTMP%NEXTBSPC)                                      40.31
            IF ( .NOT.LBFILS ) THEN                                       40.31
               FBNDFIL = BFLTMP                                           40.31
               CUBFL => FBNDFIL                                           40.31
               LBFILS = .TRUE.                                            40.31
            ELSE                                                          40.31
               CUBFL%NEXTBSPC => BFLTMP                                   40.31
               CUBFL => BFLTMP                                            40.31
            END IF                                                        40.31
            CALL ININTG('SEQ', NFSEQ, 'STA', 1)
            NBSPSS = NBSPSS + NFSEQ
          ENDIF
          DO IKO = 1, KOUNTR
            IX = CURR%JX                                                  40.31
            IF (OPTG.NE.5) IY = CURR%JY                                   40.80 40.31
            CURR => CURR%NEXTXY                                           40.31
            ALLOCATE(BGPTMP)                                              40.31
            IF (OPTG.NE.5) THEN                                           40.80
               BGPTMP%BGP(1) = KGRPNT(IX,IY)                              40.31
            ELSE                                                          40.80
               BGPTMP%BGP(1) = IX                                         40.80
            ENDIF                                                         40.80
            BGPTMP%BGP(2) = 1                                             40.31
            BGPTMP%BGP(3) = 1000                                          40.31
            BGPTMP%BGP(4) = NBSPSS                                        40.31
            BGPTMP%BGP(5) = 0                                             40.31
            BGPTMP%BGP(6) = 1                                             40.31
            NULLIFY(BGPTMP%NEXTBGP)                                       40.31
            IF ( .NOT.LBGP ) THEN                                         40.31
               FBGP = BGPTMP                                              40.31
               CUBGP => FBGP                                              40.31
               LBGP = .TRUE.                                              40.31
            ELSE                                                          40.31
               CUBGP%NEXTBGP => BGPTMP                                    40.31
               CUBGP => BGPTMP                                            40.31
            END IF                                                        40.31
          ENDDO
          NBGRPT = NBGRPT + KOUNTR
        ELSE IF (KEYWIS('VAR')) THEN
          CALL INKEYW('STA', 'PAR')
          IF (KEYWIS('PAR')) THEN
            BPARF = .TRUE.
          ELSE IF (KEYWIS('FILE')) THEN
            BPARF = .FALSE.
          ENDIF
          RLEN1 = -1.E20
          IKO = 1
          RDIST = 0.
          IBSPC1 = 1
          LFRST1 = .TRUE.
          LFRST2 = .TRUE.
          DO
            IF (LFRST1) THEN
              CALL INREAL('LEN', RLEN2, 'REQ', 0.)
              LFRST1 = .FALSE.
            ELSE
              CALL INREAL('LEN', RLEN2, 'STA', 1.E20)
            ENDIF
            IF (RLEN2.LT.0.9E20) THEN
              IF (IKO.GT.KOUNTR) THEN
                CALL MSGERR(1,
     &          'Length of segment short, boundary values ignored')       40.00
                WRITE (PRINTF, 332) RDIST, RLEN2
 332            FORMAT (' segment length=', F9.2, '; [len]=', F9.2)
              ENDIF
              IF (BPARF) THEN
                CALL INREAL ('HS',  SPPARM(1), 'REQ', 0.)
                CALL INKEYW ('STA', ' ')
                CALL INREAL ('PER', SPPARM(2), 'REQ', 0.)
                CALL INREAL ('DIR', SPPARM(3), 'REQ', 0.)
                IF (DSHAPE.EQ.1) THEN
                  CALL INREAL ('DD',  SPPARM(4), 'STA', 30.)
                  IF ((SPPARM(4).GT.360. .OR. SPPARM(4).LT. 0.).AND.      30.81
     &                .NOT.(EQREAL(SPPARM(4),OVEXCV(16)))) THEN           30.82
                    CALL MSGERR (2,'Directional spreading is less ' //    30.81
     &                             'than 0 or larger than 360 '//         30.81
     &                             'degrees and no exception value')      30.81
                  END IF                                                  30.81
                ELSE
                  CALL INREAL ('DD',  SPPARM(4), 'STA', 2.)
                  IF (SPPARM(4).LE. 0.) THEN                              30.81
                    CALL MSGERR (2,'Power of cosine is less or equal '//  30.81
     &                             'to zero')                             30.81
                  END IF                                                  30.81
                ENDIF
                NBSPEC = NBSPEC + 1
                IBSPC2 = NBSPEC
                ALLOCATE(BSTMP)                                           40.31
                BSTMP%NBS    = NBSPEC                                     40.31
                BSTMP%FSHAPE = FSHAPE                                     40.31
                BSTMP%DSHAPE = DSHAPE                                     40.31
                BSTMP%SPPARM(1:4) = SPPARM(1:4)                           40.31
                NULLIFY(BSTMP%NEXTBS)                                     40.31
                IF ( .NOT.LBS ) THEN                                      40.31
                   FBS = BSTMP                                            40.31
                   CUBS => FBS                                            40.31
                   LBS = .TRUE.                                           40.31
                ELSE                                                      40.31
                   CUBS%NEXTBS => BSTMP                                   40.31
                   CUBS => BSTMP                                          40.31
                END IF                                                    40.31
              ELSE
                IF (LFRST2) THEN
                  CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
                  LFRST2 = .FALSE.
                ELSE
                  CALL INCSTR ('FNAME', FILENM, 'STA', ' ')
                ENDIF
                IF (FILENM.NE.'    ') THEN
!                 generate new set of file data
                  NBFILS = NBFILS + 1
                  NBSPSS = NBSPEC
                  ALLOCATE(BFLTMP)                                        40.31
                  CALL BCFILE (FILENM, 'PNTS', BFLTMP, LBGP,              40.31
     &                   XCGRID, YCGRID, KGRPNT, XYTST, KGRBND,           40.31
     &                   DONALL)                                          40.05
                  IF (STPNOW()) RETURN                                    34.01
                  NULLIFY(BFLTMP%NEXTBSPC)                                40.31
                  IF ( .NOT.LBFILS ) THEN                                 40.31
                     FBNDFIL = BFLTMP                                     40.31
                     CUBFL => FBNDFIL                                     40.31
                     LBFILS = .TRUE.                                      40.31
                  ELSE                                                    40.31
                     CUBFL%NEXTBSPC => BFLTMP                             40.31
                     CUBFL => BFLTMP                                      40.31
                  END IF                                                  40.31
                ENDIF
                CALL ININTG ('SEQ', NFSEQ, 'STA', 1)
                IBSPC2 = NBSPSS + NFSEQ
                IF (IBSPC2.GT.NBSPEC)
     &                CALL MSGERR (1,'too large value for SEQ')
              ENDIF
            ELSE
              IF (IKO.GT.KOUNTR) GOTO 360                                 40.00
            ENDIF
            LFRST3 = .TRUE.
            DO
              IX = CURR%JX                                                40.31
              IF (OPTG.NE.5) THEN                                         40.80
                 IY = CURR%JY                                             40.31
                 XC2 = XCGRID(IX,IY)
                 YC2 = YCGRID(IX,IY)
              ELSE                                                        40.80
                 XC2 = xcugrd(IX)                                         40.80
                 YC2 = ycugrd(IX)                                         40.80
              ENDIF                                                       40.80
              IF (.NOT.LFRST3) THEN
                RDIST = RDIST + SQRT ((XC2-XC1)**2 + (YC2-YC1)**2)
              ENDIF
              LFRST3 = .FALSE.
              XC1 = XC2
              YC1 = YC2
              IF (RDIST.GT.RLEN2) GOTO 340
              ALLOCATE(BGPTMP)                                            40.31
              IF (OPTG.NE.5) THEN                                         40.80
                 BGPTMP%BGP(1) = KGRPNT(IX,IY)                            40.31
              ELSE                                                        40.80
                 BGPTMP%BGP(1) = IX                                       40.80
              ENDIF                                                       40.80
              BGPTMP%BGP(2) = 1                                           40.31
              W1 = (RLEN2-RDIST)/(RLEN2-RLEN1)
              BGPTMP%BGP(3) = NINT(1000.*W1)                              40.31
              BGPTMP%BGP(4) = IBSPC1                                      40.31
              BGPTMP%BGP(5) = NINT(1000.*(1.-W1))                         40.31
              BGPTMP%BGP(6) = IBSPC2                                      40.31
              NULLIFY(BGPTMP%NEXTBGP)                                     40.31
              IF ( .NOT.LBGP ) THEN                                       40.31
                 FBGP = BGPTMP                                            40.31
                 CUBGP => FBGP                                            40.31
                 LBGP = .TRUE.                                            40.31
              ELSE                                                        40.31
                 CUBGP%NEXTBGP => BGPTMP                                  40.31
                 CUBGP => BGPTMP                                          40.31
              END IF                                                      40.31
              IKO = IKO + 1
              IF (IKO.GT.KOUNTR) GOTO 340                                 40.00
              IF (.NOT.ASSOCIATED(CURR%NEXTXY)) EXIT                      40.31
              CURR => CURR%NEXTXY                                         40.31
            ENDDO
!           boundary values have been assigned, read new parameters
 340        IF (RLEN2.GT.0.9E20) GOTO 360                                 40.00
            RLEN1  = RLEN2
            IBSPC1 = IBSPC2
          ENDDO
!         update NBGRPT = number of boundary grid points
 360      NBGRPT = NBGRPT + KOUNTR                                        40.00
        ELSE
          CALL WRNKEY
        ENDIF
        IF (ASSOCIATED(TMP)) DEALLOCATE(TMP)                              40.31
      ENDIF
 900  RETURN
      END
!*********************************************************************
!                                                                    *
      SUBROUTINE BCFILE (FBCNAM, BCTYPE, BSPFIL, LBGP,                    40.31
     &                   XCGRID, YCGRID, KGRPNT,                          40.31
     &                   XYTST,  KGRBND, DONALL)                          40.31 40.05
!                                                                    *
!*********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC                                                       40.31
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
!     30.73: Nico Booij
!     30.90: IJsbrand Haagsma (Equivalence version)
!     34.01: Jeroen Adema
!     40.03, 40.13: Nico Booij
!     40.05: Ekaterini E. Kriezi
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.73, Dec. 97: New subroutine
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, June 00: function EQCSTR used to compare strings
!            July 00: option LONLAT for location coordinates introduced
!     40.05, Aug. 00: replace the source text related with the grid points
!                     interpolation coef. with a new subroutine BC_POINTS
!     40.13, Jan. 01: ! is now allowed as comment sign in a boundary file
!                     checking coordinates only for nesting situation
!                     remove declarations of unused variables
!            Nov. 01: initial size of BSPAUX array enlarged
!     40.31, Nov. 03: removing POOL-mechanism, reconsideration of this
!                     subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Nov. 04: small corrections
!
!  2. Purpose
!
!     Reads file data for boundary condition
!
!  3. Method
!
!
!  4. Argument variables
!
      REAL      XCGRID(MXC,MYC), YCGRID(MXC,MYC)
!
      INTEGER   KGRPNT(MXC,MYC)                                           40.31
      INTEGER   XYTST(*), KGRBND(*)                                       40.31
!
!       FBCNAM  char  inp    filename of boundary data file
!       BCTYPE  char  inp    if value is "NEST": nesting b.c.
!       XCGRID  real  inp    x-coordinate of computational grid points
!       YCGRID  real  inp    y-coordinate of computational grid points
!       KGRPNT  int   inp    indirect addresses of grid points
!       XYTST   int   inp    ix, iy of test points
!
!     DONALL: logic arguments declare if the nesting  boundary is open or close
!             it is defined by the users
!
      LOGICAL, INTENT(INOUT)  ::  DONALL                                  40.05
!
      CHARACTER FBCNAM *(*), BCTYPE *(*)

      TYPE(BSPCDAT) :: BSPFIL                                             40.31

      LOGICAL :: LBGP                                                     40.31
!
!  5. Parameter variables
!
!
!  6. Local variables
!
      INTEGER :: ISTATF, NDSL, NDSD, IOSTAT, IERR, NBOUNC, NANG, NFRE
      INTEGER :: IBOUNC, DORDER
      INTEGER :: IENT,IOPTT
      INTEGER :: NHEDF, NHEDT, NHEDS, IFRE , IANG
      INTEGER :: NQUANT, IQUANT, IBC, II, NBGRPT_PREV,IIPT2
      REAL    :: XP, YP, XP2, YP2
      REAL    :: FREQHZ, DIRDEG, DIRRD1,DIRRAD, EXCV
      CHARACTER BTYPE *4, HEDLIN *80
!
!    NBGRPT_PREV is the prevous number of NBGRPT
!    IIPT2 counter use for the chekinf if there are grid points on nested boundary
!
!  8. Subroutines Used
!
!       Ocean Pack command reading routines
!       SWBCPT : boundary points interpolation
!
      LOGICAL STPNOW, EQCSTR                                              40.03
!
!  9. Subroutines calling
!
!       SWREAD
!
!  10. Error messages
!
!
!  11. Remarks
!
!
!       This subroutine reads the heading of the file to determine locations
!       of boundary spectra, spectral frequencies and directions etc.
!       Reading and processing of spectral energy densities is done during
!       computation by subroutine RESPEC (file Swanmain.for)
!
!       data concerning boundary files are stored in array BFILED
!       there is a subarray for each file; it contains:
!       1.  status; 0: stationary, 1: nonstat, -1: exhausted
!       2.  time of boundary values read one before last
!       3.  time of boundary values read last
!       4.  NDSL: unit ref. num. of file containing filenames
!       5.  NDSD: unit ref. num. of file containing data
!       6.  time coding option for reading time from data file
!       8.  number of locations for which spectra are in the file
!       9.  order of reading directional information
!       10. number of spectral directions of spectra on file
!       12. number of spectral frequencies
!       14. number of heading lines per file
!       15. number of heading lines per time step
!       16. number of heading lines per spectrum
!       17. =1: energy dens., =2: variance density
!       18. =1: Cartesian direction, =2: Nautical dir.                    40.00
!       19. =1: direction spread in degr, =2: Power of Cos.               40.00
!
!
!  12. Structure
!
!       -----------------------------------------------------------------
!       Open boundary condition data file
!       Read type of file from first line of file
!       Case filetype is:
!       TPAR: make filetype TPAR
!       SWAN: make filetype SWAN
!             If b.c. type is NEST
!             Then calculate data on grid points
!             put into array BGRIDP
!             -----------------------------------------------------------
!             Read spectral directions from file into array BSPDIR
!             Read spectral frequencies from file into array BSPFRQ
!       -----------------------------------------------------------------
!       Put file characteristics into array BFILED
!       -----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL         CCOORD                                            40.03

      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'BCFILE')
!
      NDSL = 0
      IIPT2 = 0                                                           40.05
!     open data file
      NDSD = 0
      IOSTAT = 0
      CALL FOR (NDSD, FILENM, 'OF', IOSTAT)
      IF (STPNOW()) RETURN                                                34.01
!
!     --- initialize array BFILED of BSPFIL                               40.31
      BSPFIL%BFILED = 0                                                   40.31
!
!     start reading from the data file
      READ (NDSD, '(A)') HEDLIN
      IF (EQCSTR(HEDLIN,'TPAR')) THEN                                     40.03
        BTYPE  = 'TPAR'
        ISTATF = 1
        IOPTT  = 1
        NBOUNC = 1
        NANG   = 0
        NFRE   = 0
        NHEDF  = 0
        NHEDT  = 0
        NHEDS  = 0
        DORDER = 0
        ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                     40.41
        ALLOCATE(BSPFIL%BSPDIR(NANG))                                     40.41
        IF (NSTATM.EQ.0) CALL MSGERR (3,
     &      'time information not allowed in stationary mode')
        NSTATM = 1
      ELSE IF (EQCSTR(HEDLIN,'SWAN')) THEN                                40.03
        NHEDF  = 0
  10    READ (NDSD, '(A)') HEDLIN
        IF (ITEST.GE.60) WRITE (PRTEST,11) HEDLIN                         40.00
!       skip heading lines starting with comment sign ($ as in command file)
        IF (HEDLIN(1:1).EQ.COMID .OR. HEDLIN(1:1).EQ.'!') GOTO 10         40.13
        IF (EQCSTR(HEDLIN,'TIME')) THEN                                   40.03
          IF (NSTATM.EQ.0) CALL MSGERR (3,
     &    'nonstationary boundary condition not allowed '//
     &    'in stationary mode')
          NSTATM = 1
          ISTATF = 1
          BTYPE = 'SWNT'
          READ (NDSD, *) IOPTT
          READ (NDSD, '(A)') HEDLIN
          IF (ITEST.GE.60) WRITE (PRTEST,11) HEDLIN
  11      FORMAT (' heading line: ', A)
          NHEDF = 2
          NHEDT = 1
        ELSE
          ISTATF = 0
          BTYPE = 'SWNS'
          NHEDT  = 0
        ENDIF
!
!       read geographical locations
!
!       read number of boundary points
        CCOORD = .TRUE.                                                   40.03
        IF (EQCSTR(HEDLIN,'LOC')) THEN                                    40.03
          IF (BCTYPE.EQ.'NEST' .AND. KSPHER.EQ.1) CALL MSGERR (3,         40.13
     &    'Boundary locations are Cartesian, while comp. is spherical')   40.03
        ELSE IF (EQCSTR(HEDLIN,'LONLAT')) THEN                            40.03
          IF (BCTYPE.EQ.'NEST' .AND. KSPHER.EQ.0) CALL MSGERR (3,         40.13
     &    'Boundary locations are spherical, while comp. is Cartesian')   40.03
        ELSE
!         set CCOORD to False to indicate that no locations are defined   40.03
          CCOORD = .FALSE.                                                40.03
        ENDIF
        IF (CCOORD) THEN
          READ (NDSD, *) NBOUNC

          DO IBOUNC = 1, NBOUNC
            IERR = 0
            CALL REFIXY (NDSD, XP, YP, IERR)
            IF (ITEST.GE.80) THEN
              WRITE (PRTEST, *) ' B. spectrum ', IBOUNC, XP+XOFFS,
     &        YP+YOFFS, IERR                                              40.03
            ENDIF
!           in case of nesting coordinates on file are used to
!           determine interpolation coefficients
!           in other cases coordinates are ignored
            IF (BCTYPE .EQ. 'NEST') THEN
              XP2 = XP
              YP2 = YP
!
!             --- interpolate the boundaries points to the grid points of
!                 the SWAN computational grid
!
              NBGRPT_PREV = NBGRPT                                        40.05
              CALL SWBCPT (  LBGP, XCGRID, YCGRID,                        40.41 40.31 40.05
     &                       KGRPNT, XYTST,  KGRBND,XP2,YP2,IBOUNC,       40.05
     &                       NBOUNC,DONALL )                              40.05
!             check if the grid points are on nested boundary.
!             if not, stop the calculation and give an error message
              IF (NBGRPT.NE.NBGRPT_PREV) THEN                             40.05
                IIPT2 = IIPT2+1                                           40.05
              ENDIF                                                       40.05
            ENDIF
          ENDDO
!
          IF ((BCTYPE .EQ. 'NEST').AND. (IIPT2.EQ.0)) CALL MSGERR (2,
     &      'no grid points on nested boundary')                          40.05
!
          NHEDF = NHEDF + 2 + NBOUNC
          IF (ITEST.GE.60) WRITE (PRTEST,16) NBOUNC
  16      FORMAT (I6, ' boundary locations')
          READ (NDSD, '(A)') HEDLIN
          IF (ITEST.GE.60) WRITE (PRTEST,11) HEDLIN
        ELSE
          IF (BCTYPE .EQ. 'NEST') THEN
            CALL MSGERR (3, 'this file is not a true nesting file')
          ENDIF
          NBOUNC = 1
        ENDIF
!
!       read spectral resolution information
!
!       number of spectral frequencies
        IF (EQCSTR(HEDLIN(2:5),'FREQ')) THEN                              40.03
          READ (NDSD, *) NFRE
          ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                   40.31
          DO IFRE = 1, NFRE
!           read frequency in Hz and convert to radians/sec
            READ (NDSD, *) FREQHZ
            BSPFIL%BSPFRQ(IFRE) = PI2 * FREQHZ                            40.31
          ENDDO
          READ (NDSD, '(A)') HEDLIN
          IF (ITEST.GE.60) WRITE (PRTEST,11) HEDLIN                       40.00
          NHEDF = NHEDF + 2 + NFRE
        ELSE
          NFRE = 0
          IF (BCTYPE.EQ.'NEST') THEN
            CALL MSGERR (3, 'file is not a true nesting file')
          ENDIF
        ENDIF
        IF (ITEST.GE.60) WRITE (PRTEST,19) NFRE                           40.00
  19    FORMAT (I6, ' boundary frequencies')
!       number of spectral directions
        IF (EQCSTR(HEDLIN(2:4),'DIR')) THEN                               40.03
          READ (NDSD, *) NANG
          ALLOCATE(BSPFIL%BSPDIR(NANG))                                   40.31
          DO IANG = 1, NANG
!           read direction in degr and convert to radians
            READ (NDSD, *) DIRDEG
            IF (EQCSTR(HEDLIN,'N')) THEN                                  40.03
              DIRDEG = 180. + DNORTH - DIRDEG
            ENDIF
            DIRRAD = DIRDEG * PI / 180.
!           reverse order if second direction is smaller than first
            IF (IANG.EQ.1) THEN
              DIRRD1 = DIRRAD
            ELSE IF (IANG.EQ.2) THEN
              IF (DIRRAD.LT.DIRRD1) THEN
                DORDER = -1
                BSPFIL%BSPDIR(NANG) = DIRRD1                              40.31
              ELSE
                DORDER = 1
              ENDIF
              DIRRD1 = DIRRAD
            ELSE
              IF (DORDER.LT.0.) THEN
                IF (DIRRAD.GT.DIRRD1) CALL MSGERR (3,
     &          'spectral directions in file not in right order')
              ELSE
                IF (DIRRAD.LT.DIRRD1) CALL MSGERR (3,
     &          'spectral directions in file not in right order')
              ENDIF
              DIRRD1 = DIRRAD
            ENDIF
            IF (DORDER.LT.0) THEN
              BSPFIL%BSPDIR(NANG+1-IANG) = DIRRAD                         40.31
            ELSE
              BSPFIL%BSPDIR(IANG) = DIRRAD                                40.31
            ENDIF
          ENDDO
          READ (NDSD, '(A)') HEDLIN
          IF (ITEST.GE.60) WRITE (PRTEST,11) HEDLIN
          NHEDF = NHEDF + 2 + NANG
          NHEDS = 1
        ELSE
          NANG   = 0
          ALLOCATE(BSPFIL%BSPDIR(NANG))                                   40.41
          NHEDS  = 1
          DORDER = 0
        ENDIF
        IF (ITEST.GE.60) WRITE (PRTEST,23) NANG                           40.00
  23    FORMAT (I6, ' boundary directions')
!
!       read quantities (name, unit, exc. value)
!
        IF (EQCSTR(HEDLIN,'QUANT')) THEN
          READ (NDSD, *) NQUANT
          IF (.NOT.((NQUANT.EQ.1 .AND. NANG.GT.0) .OR.
     &              (NQUANT.EQ.3 .AND. NANG.EQ.0))) THEN
            CALL MSGERR (2, 'incompatible data on b.c. file')
            WRITE (PRINTF, 31) NQUANT, NANG
  31        FORMAT (I3, ' quantities; ', I5, ' directions')
          ENDIF
          DO IQUANT = 1, NQUANT
            READ (NDSD, '(A)') HEDLIN
!           if first quantity is 'EnDens' divide by Rho*Grav
            IF (IQUANT.EQ.1) THEN
              IF ( EQCSTR(HEDLIN,'ENDENS')) THEN                          40.03
!               quantity on file is energy density
                BSPFIL%BFILED(17) = 1                                     40.31
              ELSE IF ( EQCSTR(HEDLIN,'VADENS')) THEN                     40.03
!               quantity on file is variance density
                BSPFIL%BFILED(17) = 2                                     40.31
              ELSE
                CALL MSGERR (2,
     &          'Incorrect quantity in b.c.file: ' // HEDLIN(1:10))       40.03
                BSPFIL%BFILED(17) = 2                                     40.31
              ENDIF
            ELSE IF (IQUANT.EQ.2) THEN                                    40.00
!             if second quantity is 'NDIR' transform from Nautical to Cartesian dir.
              IF ( EQCSTR(HEDLIN,'NDIR')) THEN                            40.03
!               quantity on file is Nautical direction
                BSPFIL%BFILED(18) = 2                                     40.31
              ELSE IF (EQCSTR(HEDLIN,'CDIR')) THEN                        40.03
!               quantity on file is Cartesian direction
                BSPFIL%BFILED(18) = 1                                     40.31
              ELSE
                CALL MSGERR (2,
     &          'Incorrect quantity in b.c.file: ' // HEDLIN(1:10))       40.03
                BSPFIL%BFILED(18) = 1                                     40.31
              ENDIF
            ELSE IF (IQUANT.EQ.3) THEN                                    40.00
!             if third quantity is 'DSPRP' or 'POWER' power is given,
!             otherwise calculate power from dir. spread in degrees
              IF (EQCSTR(HEDLIN,'DSPRP') .OR.
     &            EQCSTR(HEDLIN,'POWER')) THEN                            40.03
!               quantity on file is power of cos
                BSPFIL%BFILED(19) = 2                                     40.31
              ELSE IF (EQCSTR(HEDLIN,'DSPR') .OR.
     &                 EQCSTR(HEDLIN,'DEGR')) THEN                        40.03
!               quantity on file is Directional spread in degr
                BSPFIL%BFILED(19) = 1                                     40.31
              ELSE
                CALL MSGERR (2,
     &          'Incorrect quantity in b.c.file: ' // HEDLIN(1:10))       40.03
                BSPFIL%BFILED(19) = 1                                     40.31
              ENDIF
            ENDIF
!           check Unit and Exception value
            READ (NDSD, '(A)') HEDLIN
            IF (IQUANT.EQ.3 .AND. EQCSTR(HEDLIN,'DEGR')) THEN
              IF (BSPFIL%BFILED(19).NE.1) THEN                            40.31
                CALL MSGERR (2, 'incompatible options in boundary file')
                BSPFIL%BFILED(19) = 1                                     40.31
              ENDIF
            ENDIF
            IF (IQUANT.EQ.1) THEN                                         40.41
               READ (NDSD, *) EXCV
               BSPFIL%BFILED(11) = NINT(EXCV)
            ELSE
               READ (NDSD, '(A)') HEDLIN
            END IF
          ENDDO
          NHEDF = NHEDF + 2 + 3*NQUANT
        ENDIF
        IF (ITEST.GE.60) WRITE (PRTEST,28) NQUANT                         40.00
  28    FORMAT (I6, ' quantities')
      ELSE
        CALL MSGERR (3, 'unsupported boundary data file')
      ENDIF
!
      ALLOCATE(BSPFIL%BSPLOC(NBOUNC))                                     40.31
      DO IBC = 1, NBOUNC
         BSPFIL%BSPLOC(IBC) = NBSPEC + IBC                                40.31
      ENDDO
      NBSPEC = NBSPEC + NBOUNC
!
!     store file reading parameters in array BFILED
!
      BSPFIL%BFILED(1)  = ISTATF                                          40.31
      BSPFIL%BFILED(2)  = -999999999                                      40.31
      BSPFIL%BFILED(3)  = -999999999                                      40.31
      BSPFIL%BFILED(4)  = NDSL                                            40.31
      BSPFIL%BFILED(5)  = NDSD                                            40.31
      BSPFIL%BFILED(6)  = IOPTT                                           40.31
      CALL COPYCH (BTYPE, 'T', BSPFIL%BFILED(7), 1, IERR)                 40.31
      BSPFIL%BFILED(8)  = NBOUNC                                          40.31
      BSPFIL%BFILED(9)  = DORDER                                          40.31
      BSPFIL%BFILED(10) = NANG                                            40.31
      BSPFIL%BFILED(12) = NFRE                                            40.31
!     ordering of data on file
      BSPFIL%BFILED(13) = 0                                               40.31
!     number of heading lines: per file, per time, per spectrum
      BSPFIL%BFILED(14) = NHEDF                                           40.31
      BSPFIL%BFILED(15) = NHEDT                                           40.31
      BSPFIL%BFILED(16) = NHEDS                                           40.31
!
      IF (ITEST.GE.80) WRITE(PRINTF,81) NBFILS, NBSPEC,
     &      (BSPFIL%BFILED(II), II=1,16)                                  40.31
  81  FORMAT (' array BFILED: ', 2I4, 2(/,8I10))
!
      RETURN
!     end of subroutine BCFILE
      END
!*********************************************************************
!                                                                    *
      SUBROUTINE BCWAMN (FBCNAM, BCTYPE, BSPFIL, LBGP,                    40.31
     &                   XCGRID, YCGRID, KGRPNT, XYTST)                   40.31
!                                                                    *
!*********************************************************************

      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC                                                       40.31
      USE M_PARALL

      IMPLICIT NONE                                                       40.13
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
!     30.73: Nico Booij
!     30.90: IJsbrand Haagsma (Equivalence version)
!     34.01: Jeroen Adema
!     40.03: Nico Booij
!     40.13: N. Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.61: Roop Lalbeharry
!
!  1. Updates
!
!     30.73, Jan. 98: new subroutine, based on older version by Weimin Luo
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Nov. 99: THD (first spectral direction in radians) added in
!                     expression for RBSDIR (directions of boundary spectrum)
!     40.03, Aug. 00: correction WAM nest with spherical SWAN
!     40.13, May  01: order of boundary points in WAM nesting file differed
!                     from order assumed in SWAN
!     40.31, Nov. 03: removing POOL-mechanism, reconsideration of this
!                     subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Nov. 06: variables USNEW, THWNEW no longer written in WAM4.5
!
!  2. PURPOSE
!
!     reads file data for WAM nesting boundary condition
!
!  3. METHOD
!
!
!  4. Argument variables
!
!       FBCNAM  char  inp    filename of boundary data file
!       BCTYPE  char  inp    if value is "NEST": nesting b.c.
!       XCGRID  real  inp    x-coordinate of computational grid points
!       YCGRID  real  inp    y-coordinate of computational grid points
!       KGRPNT  int   inp    indirect addresses of grid points
!       XYTST   int   inp    ix, iy of test points
!
!  7. Common blocks used
!
!
!  5. SUBROUTINES CALLING
!
!       SWBOUN
!
!  6. SUBROUTINES USED
!
!       Ocean Pack command reading routines
!
      LOGICAL :: STPNOW                                                   34.01


!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       data concerning boundary files are stored in array BFILED
!       there is a subarray for each file; it contains:
!       1.  status; 0: stationary, 1: nonstat, -1: exhausted
!       2.  time of boundary values read one before last
!       3.  time of boundary values read last
!       4.  NDSL: unit ref. num. of file containing filenames
!       5.  NDSD: unit ref. num. of file containing data
!       6.  time coding option for reading time from data file
!       8.  number of locations for which spectra are in the file
!       9.  order of reading directional information
!       10. number of spectral directions of spectra on file
!       12. number of spectral frequencies
!       14. number of heading lines per file
!       15. number of heading lines per time step
!       16. number of heading lines per spectrum
!       17. =1: energy dens., =2: variance density
!
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       Open file containing filnames
!       Read data file name
!       Open boundary condition data file
!       Read number of b.points, frequencies and directions
!       Generate spectral directions from file into array BSPDIR
!       Generate spectral frequencies from file into array BSPFRQ
!       For all boundary spectra do
!           read location from data file
!           transform into local cartesian or spherical coordinates       40.13
!       -----------------------------------------------------------------
!       Determine spatial step size in WAM nesting file                   40.13
!       For all spatial points in WAM file do                             40.13
!           For all other spatial points in WAM file do                   40.13
!               If the two points are neighbours                          40.13
!               Then For all computational grid points on boundary do
!                        if point is located between nest file grid points
!                        calculate interpolation coefficients
!                        and put these into array BGRIDP
!       -----------------------------------------------------------------
!       Write file characteristics into array BFILED
!       -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      INTEGER   KGRPNT(MXC,MYC), XYTST(*)
      REAL      XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      TYPE(BSPCDAT) :: BSPFIL                                             40.31
      LOGICAL :: LBGP                                                     40.31
      CHARACTER FBCNAM *(*), BCTYPE *(*)
!
!     local variables
!
      INTEGER   ISTATF, NDSL, NDSD, IOSTAT, IERR, NBOUNC, NANG, NFRE,
     &          IBOUNC, IX1, IY1, IX2, IY2, IXP, IYP, IP, MIP, INDXGR,
     &          DORDER, IOPTT
!     ISTATF    if >0 file contains nonstationary data
!     NDSL      unit ref num of namelist file
!     NDSD      unit ref num of data file
!     IOSTAT    io status
!     IERR      error status
!     NBOUNC    number of boundary locations
!     NANG      number of directions on file
!     NFRE      number of frequencies on file
!     IBOUNC    counter of boundary spectra
!     IX1
!     IY1
!     IX2
!     IY2
!     IXP
!     IYP
!     IP
!     MIP
!     INDXGR    counter of boundary grid points
!     DORDER    if <0 order of reading directions is reversed
!     IOPTT     time reading option

      INTEGER :: IIPT1=0, IIPT2=0
!     local and overall number of interpolated boundary grid points

      INTEGER :: NHEDF       ! number of heading lines at begin of file   40.13
      INTEGER :: NHEDT       ! number of heading lines per time step      40.13
      INTEGER :: NHEDS       ! number of heading lines                    40.13
      INTEGER :: IBC, IDW, ISW, IFRE, II, ISIDE, IHD    ! counters
      INTEGER :: IBNC1, IBNC2   ! counters of nesting points
      INTEGER :: IBSP1, IBSP2   ! counters of nesting points

      REAL      XP, YP, XP1, YP1, XP2, YP2, RR, RX, RY, RL2,
     &          XANG, XFRE, THD, FR1, CO, XBOU, XDELC,
     &          XLON, XLAT, XDATE, EMEAN, THQ, FMEAN, USNEW, THWNEW
!     XP        problem coordinate of a comp. grid point on the boundary
!     YP        problem coordinate of a comp. grid point on the boundary
!     XP1       problem coordinate of a boundary location
!     YP1       problem coordinate of a boundary location
!     XP2       problem coordinate of a boundary location
!     YP2       problem coordinate of a boundary location
!     RR
!     RX        vector connecting two boundary locations
!     RY        vector connecting two boundary locations
!     RL2       length **2 of vector connecting two boundary locations

      DOUBLE PRECISION, ALLOCATABLE :: XPWAM(:), YPWAM(:)
      REAL, ALLOCATABLE :: SPAUX(:)                                       40.31
      ! locations of nesting points                                       40.13
      DOUBLE PRECISION :: DXWAM, DYWAM   ! spatial step sizes in nesting file   40.13
      DOUBLE PRECISION :: DXTEST, DYTEST ! distance between two nesting points  40.13
      REAL :: DISXY          ! dim.less distance                          40.13
      REAL :: PHI            ! direction of vector (RX,RY)                40.13
      REAL :: DPHI           ! difference in direction                    40.13
      REAL :: EPS            ! tolerance                                  40.13
      REAL :: W2             ! interpolation coefficient                  40.13


      CHARACTER (LEN=4)  :: BTYPE     ! type of boundary cond.
      CHARACTER (LEN=14) :: CDATE     ! date-time
      CHARACTER (LEN=80) :: HEDLIN    ! heading line

      DOUBLE PRECISION :: DDATE, XLON0, XLAT0
!     DDATE     date-time
!     XLON0     longitude of origin of computational grid
!     XLAT0     latitude of origin of computational grid

      TYPE(BGPDAT), POINTER :: BGPTMP                                     40.31

!     subroutines used

      LOGICAL :: KEYWIS

      INTEGER, SAVE :: IENT = 0                                           40.13
      CALL STRACE (IENT, 'BCWAMN')
!
      ISTATF = 1
      IOPTT = 6
      NDSL = 0
!     open file with list of names
      CALL FOR (NDSL, FILENM,'OF',IOSTAT)
      IF (STPNOW()) RETURN                                                34.01
      READ (NDSL,'(A36)') FILENM
      CALL INKEYW ('REQ', ' ')
      IF (KEYWIS('FRE')) THEN
        BTYPE = 'WAMF'
      ELSE IF (KEYWIS('UNF')) THEN
        CALL INKEYW ('REQ', ' ')
        IF (KEYWIS('WK')) THEN
          BTYPE = 'WAMW'
        ELSE
          CALL IGNORE ('CRAY')
          BTYPE = 'WAMC'
        ENDIF
      ENDIF
!     open WAM data file
      NDSD=0
      IOSTAT = 0
      IF (BTYPE.EQ.'WAMF') THEN
        CALL FOR(NDSD,FILENM,'OF',IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
      ELSE
        CALL FOR(NDSD,FILENM,'OU',IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
      ENDIF
!
!     --- initialize array BFILED of BSPFIL                               40.31
      BSPFIL%BFILED = 0                                                   40.31
!
!     read spherical coordinates of point corresponding to
!     [xpc], [ypc] in Cartesian coordinates
!     not necessary if SWAN uses spherical coordinates                    40.03
!     if not given first point in data file is assumed
      CALL INDBLE('XGC',XLON0,'STA',-999.D0)                              40.01
      CALL INDBLE('YGC',XLAT0,'STA',-999.D0)                              40.01
      CALL ININTG('LWDATE',LWDATE,'STA',12)                               41.13
!
!     start reading from the data file
!
!     read resolution information from WAM input
!
      IF (BTYPE.EQ.'WAMF') THEN
        READ (NDSD,*) XANG, XFRE, THD, FR1, CO, XBOU, XDELC
      ELSE
!       Cray and workstation version
        READ (NDSD) XANG, XFRE, THD, FR1, CO, XBOU, XDELC
      ENDIF
!
!     number of WAM boundary points
      NBOUNC  = NINT(XBOU)
!     number of direction of WAM spectrum
      NANG = NINT(XANG)
      DORDER = -1
!     number of frequencies of WAM spectrum
      NFRE = NINT(XFRE)
!     number of heading lines: per file, per time, per spectrum
      NHEDF = 1
      NHEDT = 0
      NHEDS = 1
      IF (ITEST.GE.80) THEN
        WRITE(PRINTF,*) ' Number of frequencies in WAM:',NFRE
        WRITE(PRINTF,*) ' Number of directions in WAM:',NANG
        WRITE(PRINTF,*) ' Lowest frequency in WAM:',FR1
        WRITE(PRINTF,*) ' fi/fi-1 in WAM:',CO
      ENDIF
!
!     convert WAM wave directions to SWAN convention
!     apparently WAM uses direction TO which waves propagate !!
!
      ALLOCATE(BSPFIL%BSPDIR(NANG))                                       40.31
      DO  IDW = NANG,1,-1
        BSPFIL%BSPDIR(NANG-IDW+1) = DNORTH*DEGRAD -                       40.31 30.90
     &                 THD - REAL(IDW-1)*PI2/REAL(NANG)                   40.03
      ENDDO
      IF (ITEST.GE.50) WRITE (PRTEST,132) NANG,
     &          (BSPFIL%BSPDIR(IDW)*180./PI, IDW=1,NANG)                  40.31 30.90
 132  FORMAT (' WAMNEST dirs ', I3, (/, 20F6.0))
!
!     calculate WAM angular frequency array
!
      ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                       40.31
      BSPFIL%BSPFRQ(1) = PI2*FR1                                          40.31 30.90
      DO  ISW = 2, NFRE
        BSPFIL%BSPFRQ(ISW) = CO * BSPFIL%BSPFRQ(ISW-1)                    40.31 30.90
      ENDDO
      IF (ITEST.GE.50) WRITE (PRTEST,133) NFRE,
     &          (BSPFIL%BSPFRQ(ISW)*180./PI, ISW=1,NFRE)                  40.31 30.90
 133  FORMAT (' WAMNEST freqs ', I3, (/, 20F6.2))
      IF (NBOUNC.EQ.1) CALL MSGERR (3,
     &     'WAM nest does not work with only one nesting point')          40.13

!     allocate arrays XPWAM and YPWAM

      ALLOCATE (XPWAM(1:NBOUNC), YPWAM(1:NBOUNC))                         40.13
      ALLOCATE (SPAUX(NANG*NFRE))                                         40.31

!     read geographical locations and determine DXWAM and DYWAM           40.13

      IIPT2 = 0                                                           40.03
      DXWAM = 180.
      DYWAM = 180.
      DO IBOUNC = 1, NBOUNC
        IF (BTYPE.EQ.'WAMF') THEN
!         read boundary point coordinates from file
          READ(NDSD,*) XLON, XLAT, DDATE, EMEAN,
     &               THQ, FMEAN, USNEW, THWNEW
          IF (IBOUNC.EQ.1 .AND. ITEST.GE.80) WRITE (PRTEST, *)            40.13
     &          ' WAMNEST starting time ', DDATE                          40.13
!         read spectral densities but ignore them for the moment
          DO IFRE=1,NFRE
            READ(NDSD,*) (SPAUX(II), II=1,NANG)                           40.31 30.90
          ENDDO
        ELSE IF (BTYPE.EQ.'WAMC') THEN
!         read boundary point coordinates from file
          READ(NDSD) XLON, XLAT, XDATE, EMEAN,
     &               THQ, FMEAN, USNEW, THWNEW
          IF (IBOUNC.EQ.1 .AND. ITEST.GE.80) WRITE (PRTEST, *)            40.13
     &          ' WAMNEST starting time ', XDATE                          40.13
!         read spectral densities but ignore them for the moment
          READ(NDSD) (SPAUX(II), II=1,NANG*NFRE)                          40.31 30.90
        ELSE
!         read boundary point coordinates from file
          READ(NDSD) XLON, XLAT, CDATE(1:LWDATE), EMEAN, THQ, FMEAN       41.13
          IF (IBOUNC.EQ.1 .AND. ITEST.GE.80) WRITE (PRTEST, *)            40.13
     &          ' WAMNEST starting time ', CDATE(1:LWDATE)                41.13 40.13
!         read spectral densities but ignore them for the moment
          READ(NDSD) (SPAUX(II), II=1,NANG*NFRE)                          40.31 30.90
        ENDIF
        IF (ITEST.GE.50) WRITE (PRINTF, 178) IBOUNC, XLON, XLAT           40.13
 178    FORMAT (' boundary spectrum ', I3, ' at ', 2F12.4)                40.13
        XPWAM(IBOUNC) = XLON                                              40.13
        YPWAM(IBOUNC) = XLAT                                              40.13
!       determine DXWAM and DYWAM                                         40.13
        IF (IBOUNC.GT.1) THEN                                             40.13
          IF (ABS(XPWAM(IBOUNC)-XPWAM(IBOUNC-1)).GT.1.E-6)                40.13
     &    DXWAM = MIN (DXWAM, ABS(XPWAM(IBOUNC)-XPWAM(IBOUNC-1)))         40.13
          IF (ABS(YPWAM(IBOUNC)-YPWAM(IBOUNC-1)).GT.1.E-6)                40.13
     &    DYWAM = MIN (DYWAM, ABS(YPWAM(IBOUNC)-YPWAM(IBOUNC-1)))         40.13
        ENDIF
        IF (KSPHER.EQ.0) THEN                                             33.09
!         determine lower left corner of WAM nesting grid if not given by the user
          IF (IBOUNC.EQ.1) THEN
            IF (XLON0.LT.-900.) THEN
              XLON0 = XLON
              XLAT0 = XLAT
            ENDIF
          ENDIF
        ENDIF
      ENDDO                                                               40.13
      IF (ITEST.GE.50) WRITE (PRINTF, 182) DXWAM, DYWAM                   40.13
 182  FORMAT (' WAM step sizes: ', 2F12.4)                                40.13
      EPS = 0.01 * MIN(DXWAM,DYWAM)                                       40.13

!     determine interpolation coefficients for all couples of             40.13
!     neighbouring WAM nest points                                        40.13

      DO IBNC1 = 1, NBOUNC
 160    DO IBNC2 = IBNC1+1, NBOUNC                                        40.13
          DXTEST = ABS(XPWAM(IBNC1)-XPWAM(IBNC2))                         40.13
          DYTEST = ABS(YPWAM(IBNC1)-YPWAM(IBNC2))                         40.13
          IF ((DXTEST.LT.EPS .AND. ABS(DYTEST-DYWAM).LT.EPS) .OR.         40.13
     &        (DYTEST.LT.EPS .AND. ABS(DXTEST-DXWAM).LT.EPS)) THEN        40.13
!           points IBNC1 and IBNC2 are neighbours                         40.13
            IF (KSPHER.EQ.0) THEN                                         33.09
!             transform to local Cartesian coordinates
              XP1 = XPC + LENDEG*COS(PI*XLAT0/180.)*(XPWAM(IBNC1)-XLON0)  40.13
              YP1 = YPC + LENDEG*(YPWAM(IBNC1)-XLAT0)                     40.13
              XP2 = XPC + LENDEG*COS(PI*XLAT0/180.)*(XPWAM(IBNC2)-XLON0)  40.13
              YP2 = YPC + LENDEG*(YPWAM(IBNC2)-XLAT0)                     40.13
            ELSE
              XP1 = XPWAM(IBNC1) - XOFFS                                  40.13
              YP1 = YPWAM(IBNC1) - YOFFS                                  40.13
              XP2 = XPWAM(IBNC2) - XOFFS                                  40.13
              YP2 = YPWAM(IBNC2) - YOFFS                                  40.13
            ENDIF
!           Determine interpolation coefficients
            IBSP1 = NBSPEC+IBNC1                                          40.13
            IBSP2 = NBSPEC+IBNC2                                          40.13
            IIPT1 = 0                                                     40.03
            RX  = XP2 - XP1
            RY  = YP2 - YP1
            RL2 = RX**2 + RY**2
            IF (RL2.GT.0.) THEN
              RX  = RX/RL2
              RY  = RY/RL2
!             check whether direction of (RX,RY) corresponds to ALPC + k * 90 degr
              PHI = ATAN2(RY,RX)
              DPHI = MOD(PHI-ALPC+1.25*PI,0.5*PI)-0.25*PI
              IF (ABS(DPHI) .LT. 0.1) THEN
!               loop over boundary of comp. grid, select points between
!               (XP1,YP1) and (XP2,YP2)
                DO ISIDE = 1, 4
                  IF (ISIDE.EQ.1) THEN
                    IX1 = 1
                    IY1 = 1
                    IX2 = MXC
                    IY2 = 1
                    MIP = MXC
                  ELSE IF (ISIDE.EQ.2) THEN
                    IX1 = MXC
                    IY1 = 1
                    IX2 = MXC
                    IY2 = MYC
                    MIP = MYC
                  ELSE IF (ISIDE.EQ.3) THEN
                    IX1 = MXC
                    IY1 = MYC
                    IX2 = 1
                    IY2 = MYC
                    MIP = MXC
                  ELSE IF (ISIDE.EQ.4) THEN
                    IX1 = 1
                    IY1 = MYC
                    IX2 = 1
                    IY2 = 1
                    MIP = MYC
                  ENDIF
                  DO IP = 1, MIP-1
                    RR  = REAL(IP-1) / REAL(MIP-1)
                    IXP = IX1 + NINT(RR*REAL(IX2-IX1))
                    IYP = IY1 + NINT(RR*REAL(IY2-IY1))
                    INDXGR = KGRPNT(IXP,IYP)
                    IF (INDXGR.GT.1) THEN
                      XP = XCGRID(IXP,IYP)
                      YP = YCGRID(IXP,IYP)
!                     DISXY is relative distance from (XP,YP) to line
!                     (XP1,YP1) to (XP2,YP2)
                      DISXY = ABS(RX*(YP-YP1)-RY*(XP-XP1))
                      IF (DISXY.LT.0.1) THEN
!                       W2 is relative length of projection on line
!                       (XP1,YP1) to (XP2,YP2)
                        W2 = RX*(XP-XP1)+RY*(YP-YP1)
                        IF (W2.GT.-0.001 .AND. W2.LT.1.001) THEN
                          IF (W2.LT.0.01) W2 = 0.
                          IF (W2.GT.0.99) W2 = 1.
                          IF (ITEST.GE.80) WRITE (PRTEST, 223) IXP, IYP,
     &                          XP+XOFFS, YP+YOFFS, W2, IBNC1, IBNC2      40.13
 223                      FORMAT (' B.pnt', 2I5, 2F14.4, F6.3,
     &                            ' from ', 2I3)                          40.13
                          NBGRPT = NBGRPT + 1
                          IIPT1 = IIPT1 + 1                               40.03
                          IIPT2 = IIPT2 + 1                               40.03
                          ALLOCATE(BGPTMP)                                40.31
                          BGPTMP%BGP(1) = INDXGR                          40.31
!                         next item indicates type of boundary condition
                          BGPTMP%BGP(2) = 1                               40.31
                          BGPTMP%BGP(3) = NINT(1000. * W2)                40.31
                          BGPTMP%BGP(4) = IBSP2                           40.31
                          BGPTMP%BGP(5) = NINT(1000. * (1.-W2))           40.31
                          BGPTMP%BGP(6) = IBSP1                           40.31
                          NULLIFY(BGPTMP%NEXTBGP)                         40.31
                          IF ( .NOT.LBGP ) THEN                           40.31
                             FBGP = BGPTMP                                40.31
                             CUBGP => FBGP                                40.31
                             LBGP = .TRUE.                                40.31
                          ELSE                                            40.31
                             CUBGP%NEXTBGP => BGPTMP                      40.31
                             CUBGP => BGPTMP                              40.31
                          END IF                                          40.31
!                         test output if point is a test point
                          IF (NPTST.GT.0) THEN
                            DO IPTST = 1, NPTST
                              IF (IXP.EQ.XYTST(2*IPTST-1)+MXF-1 .AND.
     &                            IYP.EQ.XYTST(2*IPTST  )+MYF-1)
     &                        WRITE (PRTEST, 223) IXP, IYP,
     &                        XP+XOFFS, YP+YOFFS, W2, IBSP2, IBSP1
                            ENDDO
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
            IF (IIPT1.EQ.0) THEN
              WRITE (PRINTF, 218) XP1+XOFFS, YP1+YOFFS,
     &                            XP2+XOFFS, YP2+YOFFS
 218          FORMAT (' Warning: no grid points on interval from ',         40.03
     &              2F14.4, ' to ', 2F14.4)
            ENDIF
          ENDIF
        ENDDO                                                             40.13
!       first nesting point may have two valid neighbours                 40.13
        IF (IBNC1.EQ.1 .AND. IBNC2.EQ.2) GOTO 160                         40.13
      ENDDO
      IF (IIPT2.EQ.0) CALL MSGERR (2,
     &  'no grid points on nested boundary')                              40.03
      IF (ITEST.GE.60) WRITE (PRTEST,16) NBOUNC
  16  FORMAT (I6, ' boundary locations')

!     deallocate arrays XPWAM, YPWAM and SPAUX

      DEALLOCATE (XPWAM, YPWAM, SPAUX)                                    40.31 40.13

      ALLOCATE(BSPFIL%BSPLOC(NBOUNC))                                     40.31
      DO IBC = 1, NBOUNC
        BSPFIL%BSPLOC(IBC) = NBSPEC + IBC                                 40.31
      ENDDO
      NBSPEC = NBSPEC + NBOUNC
!
!     store file reading parameters in array BFILED
!
      BSPFIL%BFILED(1)  = ISTATF                                          40.31
      BSPFIL%BFILED(2)  = -999999999                                      40.31
      BSPFIL%BFILED(3)  = -999999999                                      40.31
      BSPFIL%BFILED(4)  = NDSL                                            40.31
      BSPFIL%BFILED(5)  = NDSD                                            40.31
      BSPFIL%BFILED(6)  = IOPTT                                           40.31
      CALL COPYCH (BTYPE, 'T', BSPFIL%BFILED(7), 1, IERR)                 40.31
      BSPFIL%BFILED(8)  = NBOUNC                                          40.31
      BSPFIL%BFILED(9)  = DORDER                                          40.31
      BSPFIL%BFILED(10) = NANG                                            40.31
      BSPFIL%BFILED(11) = 0                                               40.31
      BSPFIL%BFILED(12) = NFRE                                            40.31
!     ordering of data on file
      BSPFIL%BFILED(13) = 0                                               40.31
!     number of heading lines: per file, per time, per spectrum
      BSPFIL%BFILED(14) = NHEDF                                           40.31
      BSPFIL%BFILED(15) = NHEDT                                           40.31
      BSPFIL%BFILED(16) = NHEDS                                           40.31
!     quantity on file is variance density
      BSPFIL%BFILED(17) = 2                                               40.31
!
      IF (ITEST.GE.80) WRITE(PRINTF,81) NBFILS, NBSPEC,
     &      (BSPFIL%BFILED(II), II=1,16)                                  40.31
  81  FORMAT (' array BFILED: ', 2I4, 2(/,8I10))
!
!     Rewind input file for proper start
      REWIND (NDSD)
!     read heading line
      IF (BTYPE.EQ.'WAMF') THEN
        DO IHD = 1, BSPFIL%BFILED(14)                                     40.31
          READ (NDSD, '(A)') HEDLIN
          IF (ITEST.GE.80) WRITE (PRINTF, 212) HEDLIN
 212      FORMAT (' heading line: ', A)
        ENDDO
      ELSE
        DO IHD = 1, BSPFIL%BFILED(14)                                     40.31
          READ (NDSD)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE BCWAMN

!*********************************************************************
!                                                                    *
      SUBROUTINE BCWW3N (FBCNAM, BCTYPE, BSPFIL, LBGP,                    40.31
     &                   XCGRID, YCGRID, KGRPNT,                          40.31
     &                   XYTST,  KGRBND, DONALL)                          40.31
!                                                                    *
!*********************************************************************
!
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC                                                       40.31
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
!     40.05 : Ekaterini E. Kriezi
!     40.13 : Nico Booij
!     40.31 : Marcel Zijlema
!     40.41 : Marcel Zijlema
!
!  1. Updates
!
!     40.05, Aug. 00: new subroutine
!     40.13, Jan. 01: remove declarations of unused variables
!     40.31, Nov. 03: removing POOL-mechanism, reconsideration this
!                     subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     reads file data for WaveWatch III boundary condition
!
!  3. Method
!
!      open boundaries files
!      read from ASSCII files the points where the energy desity is given and
!      interpolate them to the grid points of the SWAN computational grid
!
!  4. Argument variables
!
      INTEGER, INTENT(IN)     ::  KGRPNT(MXC,MYC)
!                                 indirect addresses of computational grid points
      INTEGER, INTENT(IN)     ::  KGRBND(*)
!                                 array of boundary grid points
      INTEGER, INTENT(IN)     ::  XYTST(*)
!                                 array of (ix,iy) of test points
      REAL, INTENT(IN)        ::  XCGRID(MXC,MYC), YCGRID(MXC,MYC)
!                                 coordinates of computational grid points
!     FBCNAM  char  inp    filename of boundary data file
!     BCTYPE  char  inp    boundary condition type, is 'WW3N' in this case
!
      CHARACTER FBCNAM *(*), BCTYPE *(*)
!
!     DONALL : logic arguments declare if the boundary is open or close
!
      LOGICAL   :: DONALL

      TYPE(BSPCDAT) :: BSPFIL                                             40.31
      LOGICAL :: LBGP                                                     40.31
!
!  5. Parameter variables
!
!     --
!
!  6. Local variables
!
!     IENT         number of entries into this subroutine
!     IBC          spectrum counter
!
      INTEGER   :: IENT, IBC
!
      REAL,    ALLOCATABLE :: FRQ_ARRAY(:), DIR_ARRAY(:)
!
!     IHD          counter of heading lines
!     WWDATE       date in boundary file
!     WWTIME       time in boundary file
!
      INTEGER   :: WWDATE, WWTIME
!
!     ISTATF    if >0 file contains nonstationary data
!     NDSL      unit ref num of namelist file
!     NDSD      unit ref num of data file
!     IOSTAT    io status
!     IERR      error status
!     NBOUNC    number of boundary locations
!     NANG      number of directions on file
!     NFRE      number of frequencies on file
!     DORDER    if <0 order of reading directions is reversed
!     IOPTT     time reading option
!     IBOUNC    counter of boundary points
!     IGRBND    counter of boundary grid(swan grid)  points
!     II        counter
!     NHEDF     number of heading lines per file
!     NHEDT     number of heading lines per time step
!     NHEDS     number of heading lines per spectrum
!
      INTEGER            :: ISTATF, NDSL, NDSD, IOSTAT, IERR
      INTEGER            :: NBOUNC, NANG, NFRE
      INTEGER            :: IBOUNC
      INTEGER            :: DORDER, IOPTT
      INTEGER            :: NHEDF, NHEDT, NHEDS, NBGRPT_PREV
      INTEGER            :: IFRE, IANG,II, IIPT2
!
!     DUM_A     real number used for reading a file but not used in any calculation
!     XLON      longitude
!     XLAT      latitude
!     XP2       problem coordinate of a boundary location
!     YP2       problem coordinate of a boundary location
!     DIRRD1
!     NBGRPT_PREV is the prevous number of NBGRPT
!     IIPT2 counter use for the chekinf if there are grid points on nested boundary
!
      REAL               :: DUM_A, XLON, XLAT,XP2,YP2,DIRRD1
!
      CHARACTER (LEN=4)  :: BTYPE
!                           type of boundary cond.
      CHARACTER (LEN=24) :: HEDLINT
!                           WW3 version
      CHARACTER (LEN=30) :: GNAME
!                           name of test case readed from b. file
      CHARACTER (LEN=12) :: PTNME
!                           name of b. point
!     XLON0     longitude of origin of computational grid
!     XLAT0     latitude of origin of computational grid
!
      DOUBLE PRECISION   :: XLON0, XLAT0
!
!       FBCNAM  char  inp    filename of boundary data file
!       BCTYPE  char  inp    if value is "NEST": nesting b.c.
!       XCGRID  real  inp    x-coordinate of computational grid points
!       YCGRID  real  inp    y-coordinate of computational grid points
!       KGRPNT  int   inp    indirect addresses of grid points
!       XYTST   int   inp    ix, iy of test points
!
!  8. Subroutines used
!
!     Ocean Pack command reading routines
!     SWBCPT, STPNOW
!
      LOGICAL   :: STPNOW
!
!  9. Subroutines calling
!
!     SWREAD
!
! 10. Error messages
!
!     ---
!
!  11. Remarks
!
!       data concerning boundary files are stored in array BFILED
!       there is a subarray for each file; it contains:
!       1.  status; 0: stationary, 1: nonstat, -1: exhausted
!       2.  time of boundary values read one before last
!       3.  time of boundary values read last
!       4.  NDSL: unit ref. num. of file containing filenames
!       5.  NDSD: unit ref. num. of file containing data
!       6.  time coding option for reading time from data file
!       8.  number of locations for which spectra are in the file
!       9.  order of reading directional information
!       10. number of spectral directions of spectra on file
!       12. number of spectral frequencies
!       14. number of heading lines per file
!       15. number of heading lines per time step
!       16. number of heading lines per spectrum
!       17. =1: energy dens., =2: variance density, =3 variance energy density (k)
!       18. =1: Cartesian direction, =2: Nautical dir.
!       19. =1: direction spread in degr, =2: Power of Cos.
!       20.  depth of boundary points
!
!  12. Structure
!
!       -----------------------------------------------------------------
!       Open boundary condition data file
!       Read type of file from first line of file
!       If the headline is WAVEWATCH III SPECTRA
!
!       then  b.c. type is WW3N
!       -----------------------------------------------------------
!          Read spectral directions from file and write them into
!          array BSPDIR
!          Read spectral frequencies from fileand write them into
!          array BSPFRQ
!          For all boubdaries points do
!            read location from data file
!            transform into local cartesian coordinates (if nesesery)
!            Then calculate data on grid points
!
!       -----------------------------------------------------------------
!       Put file characteristics into array BFILED
!       -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'BCWW3N')
!
!     NDSL unit ref number for namelist files
      NDSL = 0
      ISTATF = 1
!     number of heading lines: per file, per time, per spectrum
      NHEDF = 0
      NHEDT = 0
      NHEDS = 0
      DORDER  = -1
      IOPTT = 1
      IIPT2 = 0
!
!     open data file NDSD unit ref number for data files
      NDSD = 0
      IOSTAT = 0
!
      CALL FOR (NDSD, FILENM , 'OF', IOSTAT)
      IF (STPNOW()) RETURN
!
!     --- initialize array BFILED of BSPFIL                               40.31
      BSPFIL%BFILED = 0                                                   40.31
!
!     start reading from the boundary data file
!     HEDLINT = 'WAVEWATCH III SPECTRA' = header of the file
!     read from boundary file the header , number of frequencies NFR,
!     number of direction NANG,  number of boundaries points NBOUNC,
!     Name of the points  GNAME
!
      READ (NDSD, 1944) HEDLINT,NFRE, NANG, NBOUNC, GNAME
!
      IF (HEDLINT(2:22) .NE. 'WAVEWATCH III SPECTRA')
     &  CALL MSGERR (3, 'file is not a WW3 spectral file')
      IF (NBOUNC.LT.2)  CALL MSGERR
     &  (3, 'SWAN need at least 2 boundary points for nesting')
!
!     ISTATF related with stationary and non stationary mode
!
      BTYPE = 'WW3N'
!
!     read frequencies from WW3 boundary file
!
      ALLOCATE (FRQ_ARRAY(1:NFRE), DIR_ARRAY(1:NANG) )
!
!     read frequency
!     FRQ_ARRAY(IFRE) =   SIG(IK)/(2*PI)
!
      READ (NDSD,1945) (FRQ_ARRAY(IFRE) ,IFRE=1,NFRE)
!
      ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                       40.31
      DO IFRE = 1, NFRE
        BSPFIL%BSPFRQ(IFRE) =  FRQ_ARRAY(IFRE)*2*PI                       40.31
      ENDDO
!
      IF (ITEST.GE.60) THEN
        WRITE(PRTEST,*) ' HEDLINT ',' NFRE ',' NANG ',' NBOUNC ',
     &                  ' GNAME'
        WRITE(PRTEST,1944) HEDLINT,NFRE, NANG, NBOUNC, GNAME
        WRITE (PRTEST,*) 'Frequencies read from boundary file ', FILENM
        WRITE (PRTEST,*) (FRQ_ARRAY(IFRE),IFRE = 1,NFRE)
      ENDIF
!
!     read direction from WW3 boundary file
!     DIR_ARRAY(IANG) = MOD(2.5*PI-TH(ITH),TPI)
!     there are in radians but is not in right order related to SWAN
!
      READ (NDSD,1946) (DIR_ARRAY(IANG),IANG=1,NANG)
!
!     put values in right order. The value of the DIR_ARRAY(i) should be
!     smaller that the DIR_ARRAY(i-1)
!     in the opposite situation make DIR_ARRAY(i) = DIR_ARRAY(i) - 2*PI
!
      DIR_ARRAY(:) = PI*DNORTH/180 - DIR_ARRAY(:)                         40.15

      ALLOCATE(BSPFIL%BSPDIR(NANG))                                       40.31
      DO IANG = 1, NANG
        IF (IANG.EQ.1) THEN
          BSPFIL%BSPDIR(1) = DIR_ARRAY(IANG)                              40.31
          DIRRD1 = BSPFIL%BSPDIR(1)                                       40.31
        ELSE
          IF (DIR_ARRAY(IANG).LT.DIRRD1) THEN
            BSPFIL%BSPDIR(IANG) = 2*PI+DIR_ARRAY(IANG)                    40.31
            DIRRD1 = BSPFIL%BSPDIR(IANG)                                  40.31
          ELSE
            BSPFIL%BSPDIR(IANG) = DIR_ARRAY(IANG)                         40.31
            DIRRD1 =  BSPFIL%BSPDIR(IANG)                                 40.31
          ENDIF
        ENDIF
      ENDDO

      IF(ITEST.GE.60) THEN
        WRITE (PRTEST,*) 'Directions read from boundary file ',
     &                    FILENM
        WRITE (PRTEST,1946) (DIR_ARRAY(IANG),IANG = 1,NANG)
      ENDIF
!
!     Time
      READ (NDSD, 900) WWDATE,WWTIME
!
!     Read from boundary file info about the boundary points(b.p): name of b. p.,
!     geographical location of b.p., depth, wind u-velocity  and direction at the b.p.
!     current velocity and direction at the b.p.
!
!     If  DONALL = .TRUE. boundary data correspond to an open boundary otherwise
!     it is continue the interpolation of the grid point between the last and the
!     first point
!
      DO IBOUNC = 1, NBOUNC
        IERR = 0
!
!       latitude =  XLAT
!       longitude = XLON
!       A real which is not used in the computation
!
        READ (NDSD,901) PTNME, XLAT, XLON, DUM_A, DUM_A,
     &                   DUM_A, DUM_A, DUM_A
!       Pass over the lines where the energy spectra is written in the boundary file.
!       The energy spectra is going to be read later, in the subroutine RESPEC
!
        READ (NDSD,902) ((DUM_A, IFRE = 1,NFRE),IANG = 1,NANG)
!
        IF (ITEST.GE.80) THEN
          WRITE (PRTEST, *) ' B. spectrum WW3 ', IBOUNC, XLON,
     &    XLAT, IERR
        ENDIF
!
!       in case of nesting coordinates on file are used to determine interpolation
!       coefficients
!
        IF (KSPHER.EQ.0) THEN
!
!       if SWAN uses Cartesian coordinates, then transform the spherical coordinates
!       of the boundary point to local Cartesian coordinates
!
          IF (IBOUNC.EQ.1) THEN
            CALL INDBLE('XGC',XLON0,'REQ',-999.D0)
            CALL INDBLE('YGC',XLAT0,'REQ',-999.D0)
            IF (XLON0.LT.-900.) THEN
              XLON0 = (XOFFS -XPC)/LENDEG
              XLAT0 = (YOFFS-YPC)/LENDEG
            ENDIF
          ENDIF
!
          XP2 = XPC + LENDEG*COS(PI*XLAT0/180.)*(XLON-XLON0)
          YP2 = YPC + LENDEG*(XLAT-XLAT0)
!
        ELSE
          XP2 = XLON-XOFFS
          YP2 = XLAT-YOFFS
        ENDIF
!
!       --- interpolate the boundaries points to the grid points of
!           the SWAN computational grid
!
        NBGRPT_PREV = NBGRPT
        CALL SWBCPT (  LBGP, XCGRID, YCGRID,                              40.41 40.31
     &                 KGRPNT, XYTST,  KGRBND,XP2,YP2,IBOUNC,
     &                 NBOUNC, DONALL )
!       check if the grid points are on nested boundary.
!       if not, stop the calculation and give an error message
        IF (NBGRPT.NE.NBGRPT_PREV) THEN
          IIPT2 = IIPT2+1
        ENDIF
      ENDDO
!
      IF (IIPT2.EQ.0) CALL MSGERR (2,
     &  'no grid points on nested boundary')
!
      IF (ITEST.GE.60) WRITE (PRTEST,16) NBOUNC
  16  FORMAT (I6, ' boundary locations')
!
!     quantity on file is energy density
      BSPFIL%BFILED(17) = 1                                               40.31
!
      NHEDF = NHEDF + CEILING(NFRE/8.) + CEILING(NANG/7.)+1
      NHEDS = 2
!     NHEDT: calculated in the RBFILE subroutine for each time step
!
      ALLOCATE(BSPFIL%BSPLOC(NBOUNC))                                     40.31
      DO IBC = 1, NBOUNC
        BSPFIL%BSPLOC(IBC) = NBSPEC + IBC                                 40.31
      ENDDO
!
      NBSPEC = NBSPEC + NBOUNC
!
!     store file reading parameters in array BFILED
!
      BSPFIL%BFILED(1)  = ISTATF                                          40.31
      BSPFIL%BFILED(2)  = -999999999                                      40.31
      BSPFIL%BFILED(3)  = -999999999                                      40.31
      BSPFIL%BFILED(4)  = NDSL                                            40.31
      BSPFIL%BFILED(5)  = NDSD                                            40.31
      BSPFIL%BFILED(6)  = IOPTT                                           40.31
      CALL COPYCH (BTYPE, 'T', BSPFIL%BFILED(7), 1, IERR)                 40.31
      BSPFIL%BFILED(8)  = NBOUNC                                          40.31
      BSPFIL%BFILED(9)  = DORDER                                          40.31
      BSPFIL%BFILED(10) = NANG                                            40.31
      BSPFIL%BFILED(11) = 0                                               40.31
      BSPFIL%BFILED(12) = NFRE                                            40.31
!     ordering of data on file
      BSPFIL%BFILED(13) = 0                                               40.31
!     number of heading lines: per file, per time, per spectrum
      BSPFIL%BFILED(14) = NHEDF                                           40.31
      BSPFIL%BFILED(15) = NHEDT                                           40.31
      BSPFIL%BFILED(16) = NHEDS                                           40.31
!     quantity on file is energy density (k)
      BSPFIL%BFILED(17) = 3                                               40.31
!
      IF (ITEST.GE.80) WRITE(PRINTF,81) NBFILS, NBSPEC,
     &      (BSPFIL%BFILED(II), II=1,16)                                  40.31
!
!      Rewind input file for proper start
      REWIND (NDSD)

  81  FORMAT (' array BFILED: ', 2I4, 2(/,8I10))
  900 FORMAT (I8.8,I7.6)
  901 FORMAT (A12,2F7.2,F10.1,2(F7.2,F6.1))
  902 FORMAT (7E11.3)
 1944 FORMAT (A23,1X,3I6,1X,A33)
 1945 FORMAT (8E10.3)
 1946 FORMAT (7E11.3)

      DEALLOCATE (FRQ_ARRAY,DIR_ARRAY)

      RETURN

      END SUBROUTINE BCWW3N

!***********************************************************************
!
      SUBROUTINE SWBCPT ( LBGP, XCGRID, YCGRID,                           40.41 40.31
     &                    KGRPNT, XYTST,  KGRBND,XP2,YP2,IBOUNC,
     &                    NBOUNC,DONALL )
!
!************************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_BNDSPEC
      USE M_PARALL
      USE SwanGriddata                                                    40.80
      USE SwanGridobjects                                                 40.80
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
!     30.73, 40.13: Nico Booij
!     40.05: Ekaterini E. Kriezi
!     40.31: Tim Campbell and John Cazes
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     41.14: Nico Booij
!
!  1. Updates
!
!     40.05, Aug. 00: remove the code to a separate subroutine
!     40.13, Jan. 01: remove declarations of unused variables
!     40.31, Jul. 03: initializations XP0, XP1, YP0, YP1
!     40.31, Nov. 03: removing POOL mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.14, Jul. 10: error in nesting unstructured grid corrected
!
!  2. Purpose
!
!     interpolation of own boundary points to nesting boundary points
!
!  3. Method
!
!     calculate interpolation coefficients between the nesting boundary grid and
!     the SWAN computational grid
!
!  4. Argument variables
!
!     NBOUNC:     max number of boundaries points
!     IBOUNC:     counter
!
!
      INTEGER, INTENT(IN)     ::  KGRPNT(MXC,MYC)  ! indirect addresses of computational grid points
      INTEGER, INTENT(IN)     ::  KGRBND(*)        ! array of boundary grid points
      INTEGER, INTENT(IN)     ::  XYTST(*)         ! array of (ix,iy) of test points
!
      INTEGER, INTENT(IN)     ::  IBOUNC,NBOUNC
!
      REAL                    ::  XP2, YP2         ! coordinates of a point in spectral file   41.14
      REAL, INTENT(IN)        ::  XCGRID(MXC,MYC), YCGRID(MXC,MYC)  ! coordinates of computational grid points
!
!     DONALL : logic arguments declare if the nesting boundary is open or close
!              it is defined by the users
!
      LOGICAL, INTENT(INOUT)  ::  DONALL

      LOGICAL :: LBGP                                                     40.31
!
!  5. Parameter variables
!
!     --
!
!  6. Local variables
!
!     XP1       problem coordinate of a boundary location
!     YP1       problem coordinate of a boundary location
!     XP2       problem coordinate of a boundary location
!     YP2       problem coordinate of a boundary location
!     DISXY     distance in (x,y)-space
!     W2        is relative length of projection on line
!     IIPT1     counter checking the grid points related to the nesting boundary
!
      INTEGER, SAVE      :: IBSP0 = 1, IBSP1 = 1, IENT = 0
      INTEGER            :: IBSP2,IGRBND,INDXGR
      INTEGER            :: IXP,IYP,IIPT1
!
      REAL, SAVE         :: XP0=0., XP1=0., YP0=0., YP1=0.                40.31
      REAL               :: XP, YP, RX, RY, RL2
      REAL               :: DISXY, W2
      TYPE(BGPDAT), POINTER :: BGPTMP                                     40.31
!
      TYPE(verttype), DIMENSION(:), POINTER :: vert                       40.80
!
!  7. Common blocks used
!
!
!  8. Subroutines Used
!
!
!  9. Subroutines calling
!
!      BCFILE : open and read Swan nesting files
!      BCWW3N : open and read WW3 nesting files
!
!  10. Error messages
!
!       ---
!
!  11. Remarks
!
!  12. Structure
!
!      For all computational grid points on boundary do
!          if point is located between nest file grid points
!             calculate interpolation coefficients
!          if DONALL is TRUE
!             the nesting boundary remain open
!          else DONALL is FALSE (default case)
!             boundary is close, it do interpolation between the last and the first point
!             put interpolation coefficients into array BGRIDP
!
!  13. Source text
!
      CALL STRACE (IENT, 'SWBCPT')                                        40.41
!
      IF (NBOUNC.EQ.1) CALL MSGERR (2,                                    41.14
     &  'Nesting procedure does not work if file has only 1 spectrum')    41.14
!
!     point to vertices
!
      vert => gridobject%vert_grid                                        40.80
!
      IBSP2 = NBSPEC+IBOUNC
      IIPT1 = 0
!
 201  IF (IBOUNC.EQ.1) THEN
!        first point found in a spectral input file                       41.14
         XP0   = XP2
         YP0   = YP2
         IBSP0 = IBSP2
         IIPT1 = 1                                                        40.41
      ELSE
!        RX, RY difference vector between consecutive points of spectral file  41.14
         RX  = XP2 - XP1
         RY  = YP2 - YP1
         RL2 = RX**2 + RY**2
         IF (RL2.GT.0.) THEN
           RX  = RX/RL2
           RY  = RY/RL2
!
!          loop over boundary of comp. grid, select boundary points
!          between (XP1,YP1) and (XP2,YP2)
!
           IF (OPTG.EQ.5) THEN                                            40.80
             DO IXP = 1, nverts                                           40.80
               IF ( vert(IXP)%atti(VMARKER) == 1 .AND.                    40.80
     &              vert(IXP)%atti(VBC) == 0 ) THEN                       40.80
                 XP = vert(IXP)%attr(VERTX)                               41.14
                 YP = vert(IXP)%attr(VERTY)                               41.14
!
!                DISXY is relative distance from (XP,YP) to line
!                (XP1,YP1) to (XP2,YP2) with respect to the
!                length of that line
!
                 DISXY = ABS(RX*(YP-YP1)-RY*(XP-XP1))                     41.14
!
                 IF (DISXY.LT.0.1) THEN
!                  W2 is relative length of projection on line
!                  (XP1,YP1) to (XP2,YP2)
                   W2 = RX*(XP-XP1)+RY*(YP-YP1)                           41.14
                   IF (W2.GT.-0.001 .AND. W2.LT.1.001) THEN               41.14
                     IF (W2.LT.0.01) W2 = 0.                              41.14
                     IF (W2.GT.0.99) W2 = 1.                              41.14
                     IF (ITEST.GE.80) WRITE (PRTEST, *) ' B.pnt',
     &                        IXP, XP, YP, W2, IBSP2, 1.-W2, IBSP1        41.14
                     NBGRPT = NBGRPT + 1
                     IIPT1 = IIPT1 + 1
!
                     ALLOCATE(BGPTMP)                                     40.31
                     BGPTMP%BGP(1) = IXP                                  40.31
!                    next item indicates type of boundary condition
                     BGPTMP%BGP(2) = 1                                    40.31
                     BGPTMP%BGP(3) = NINT(1000. * W2)                     40.31
                     BGPTMP%BGP(4) = IBSP2                                40.31
                     BGPTMP%BGP(5) = NINT(1000. * (1.-W2))                40.31
                     BGPTMP%BGP(6) = IBSP1                                40.31
                     vert(IXP)%atti(VBC) = 1                              40.80
                     NULLIFY(BGPTMP%NEXTBGP)                              40.80
                     IF ( .NOT.LBGP ) THEN                                40.80
                        FBGP = BGPTMP                                     40.80
                        CUBGP => FBGP                                     40.80
                        LBGP = .TRUE.                                     40.80
                     ELSE                                                 40.80
                        CUBGP%NEXTBGP => BGPTMP                           40.80
                        CUBGP => BGPTMP                                   40.80
                     ENDIF                                                40.80
                   ENDIF                                                  41.14
                 ENDIF                                                    41.14
               ENDIF                                                      40.80
             ENDDO                                                        40.80
           ELSE
!
!           KGRBND grid addresses on boundary points, NGRBND number of grid points
!           on computational grid boundary
!
             DO IGRBND = 1, NGRBND
               IXP = KGRBND(2*IGRBND-1)
               IYP = KGRBND(2*IGRBND)
!
               IF (IXP.GT.0 .AND.IYP.GT.0) THEN
!
!                (IXP,IYP) is a valid boundary point
!
                 INDXGR = KGRPNT(IXP,IYP)
                 XP = XCGRID(IXP,IYP)
                 YP = YCGRID(IXP,IYP)
!
!                DISXY is relative distance from (XP,YP) to line
!                (XP1,YP1) to (XP2,YP2) with respect to the
!                length of that line
!
                 DISXY = ABS(RX*(YP-YP1)-RY*(XP-XP1))
!
                 IF (DISXY.LT.0.1) THEN
!                  W2 is relative length of projection on line
!                  (XP1,YP1) to (XP2,YP2)
                   W2 = RX*(XP-XP1)+RY*(YP-YP1)
                   IF (W2.GT.-0.001 .AND. W2.LT.1.001) THEN
                     IF (W2.LT.0.01) W2 = 0.
                     IF (W2.GT.0.99) W2 = 1.
                     IF (ITEST.GE.80) WRITE (PRTEST, *) ' B.pnt',
     &                        IXP, IYP, XP, YP, W2, IBSP2, 1.-W2, IBSP1   40.41
                     NBGRPT = NBGRPT + 1
                     IIPT1 = IIPT1 + 1
!
                     ALLOCATE(BGPTMP)                                     40.31
                     BGPTMP%BGP(1) = INDXGR                               40.31
!                    next item indicates type of boundary condition
                     BGPTMP%BGP(2) = 1                                    40.31
                     BGPTMP%BGP(3) = NINT(1000. * W2)                     40.31
                     BGPTMP%BGP(4) = IBSP2                                40.31
                     BGPTMP%BGP(5) = NINT(1000. * (1.-W2))                40.31
                     BGPTMP%BGP(6) = IBSP1                                40.31
                     NULLIFY(BGPTMP%NEXTBGP)                              40.31
                     IF ( .NOT.LBGP ) THEN                                40.31
                        FBGP = BGPTMP                                     40.31
                        CUBGP => FBGP                                     40.31
                        LBGP = .TRUE.                                     40.31
                     ELSE                                                 40.31
                        CUBGP%NEXTBGP => BGPTMP                           40.31
                        CUBGP => BGPTMP                                   40.31
                     END IF                                               40.31
!
!                    test output if point is a test point
!
                     IF (NPTST.GT.0) THEN
                       DO IPTST = 1, NPTST
                         IF (IXP.EQ.XYTST(2*IPTST-1)+MXF-1 .AND.
     &                       IYP.EQ.XYTST(2*IPTST  )+MYF-1)
     &                      WRITE (PRTEST, 223)
     &                      IXP-1,IYP-1, XP+XOFFS, YP+YOFFS, W2, IBSP2,
     &                      IBSP1
 223                        FORMAT (' B.pnt', 2I5, 2F9.0, F6.3, 2I3)
                       ENDDO
                     ENDIF
                   ENDIF
                 ENDIF
               ENDIF
             ENDDO
           ENDIF
         ENDIF
      ENDIF
!
      IF (IIPT1.EQ.0) THEN
         WRITE (PRINTF, 218) XP1+XOFFS, YP1+YOFFS,
     &   XP2+XOFFS, YP2+YOFFS
 218     FORMAT (' Warning: no grid points on interval from ', 2F12.4,
     &           ' to ', 2F12.4)
      ENDIF
!
      XP1   = XP2
      YP1   = YP2
      IBSP1 = IBSP2
!
      IF (IBOUNC.EQ.NBOUNC) THEN
         IF (.NOT. DONALL) THEN
!
!           process grid points between last and first boundary point
!
            DONALL = .TRUE.
            XP2    = XP0
            YP2    = YP0
            IBSP2  = IBSP0
            GOTO 201
         ENDIF
      ENDIF

      RETURN
!
      END SUBROUTINE SWBCPT
!
!*********************************************************************
!                                                                    *
      LOGICAL FUNCTION BOUNPT (IX,IY,KGRPNT)
!                                                                    *
!*********************************************************************
!
      USE SWCOMM3                                                         40.41
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
!     40.00, 40.03: Nico Booij
!
!  1. UPDATE
!
!       Feb. 1998, ver. 40.00: new subroutine
!       40.03, Sep. 00: inconsistency with manual corrected
!
!  2. PURPOSE
!
!       determine whether a grid point is a point where a boundary condition
!       can be applied
!
!  3. METHOD
!
!
!  4. PARAMETERLIST
!
!       IX, IY  int   inp    grid point indices
!       KGRPNT  int   inp    indirect addresses of grid points
!
!  5. SUBROUTINES CALLING
!
!       BCFILE
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
!     KGRPNT(IX,IY)=1 means that (IX,IY) is not an active grid point
!
!  9. STRUCTURE
!
!     -----------------------------------------------------------------
!     Make BOUNPT = False
!     If the grid point is not active
!     Then return
!     -----------------------------------------------------------------
!     If grid point is on the outer boundary
!     Then make BOUNPT = True
!          return
!     -----------------------------------------------------------------
!     If a neighbouring grid point is inactive
!     Then make BOUNPT = True
!          return
!     -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      INTEGER IX,IY,KGRPNT(MXC,MYC)
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT, 'BOUNPT')
!
      BOUNPT = .FALSE.
!
      IF (IX.LE.0)   RETURN
      IF (IY.LE.0)   RETURN
      IF (IX.GT.MXC) RETURN
      IF (IY.GT.MYC) RETURN
!
!     If the grid point is not active
!     Then return
!
      IF (KGRPNT(IX,IY).LE.1) RETURN
!
!     If grid point is on the outer boundary
!     Then make BOUNPT = True
!          return
!
      IF (IX.EQ.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IX.EQ.MXC) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IY.EQ.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IY.EQ.MYC) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
!
!     If a neighbouring grid point is inactive
!     Then make BOUNPT = True
!          return
!
      IF (KGRPNT(IX-1,IY).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX+1,IY).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX,IY-1).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX,IY+1).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      RETURN
      END
!*********************************************************************
!                                                                    *
      SUBROUTINE RETSTP (LXYTST, XYTST, KGRPNT, KGRBND, XCGRID, YCGRID,   40.80
     &                   SPCSIG, SPCDIR)                                  40.31
!                                                                    *
!*********************************************************************
!
      USE OCPCOMM2                                                        40.95
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.31
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
!     30.82: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     34.01: Jeroen Adema
!     40.04: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Apr. 98: New subroutine
!     30.82, Oct. 98: Updated description of several variables
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Dec. 99: XOFFS and YOFFS introduced in write statements
!            May  00: ITMOPT written to heading of file instead of 1
!     40.04, Aug. 00: Added error message if testpoints are defined
!                     before bottom grid is read
!     40.13, Jan. 01: two output strings corrected
!            May  01: two incorrect units changed from m2/2 to m2/s
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Dec. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Jun. 07: extension to unstructured grids
!
!  2. PURPOSE
!
!     read test points, generate output point set 'TESTPNTS',
!     read source term filenames
!
!  3. METHOD
!
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
! i   LXYTST: Maximum length of array XYTST                               40.80
! i   XYTST : Grid point indices of test points                           30.82
!
!     MPTST : Maximum number of test points                               30.82
!
      INTEGER LXYTST, MPTST                                               40.80 30.82
      INTEGER XYTST(LXYTST)                                               40.80 30.82
!
!  5. SUBROUTINES CALLING
!
!     SWREAD
!
!  6. SUBROUTINES USED
!
      LOGICAL STPNOW                                                      34.01
!
!  7. Common blocks used
!
!
!  8. REMARKS
!
!
!  9. STRUCTURE
!
!     -----------------------------------------------------------------
!     Repeat
!         read (i,j) identifying a test point
!         store values in array XYTST
!     -----------------------------------------------------------------
!     Generate output point set 'TESTPNTS'
!     write coordinates of test points into array OUTDA
!     -----------------------------------------------------------------
!     If 1D output of source terms is requested
!     Then open file
!          write general data into the file
!     -----------------------------------------------------------------
!     If 2D output of source terms is requested
!     Then open file
!          write general data into the file
!     -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      INTEGER   KGRPNT(MXC,MYC), KGRBND(*)                                40.31
      LOGICAL   KEYWIS, LOCGRI                                            40.00
      TYPE(OPSDAT), POINTER :: OPSTMP                                     40.31
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'RETSTP')
!
      IF (OPTG.NE.5) THEN                                                 40.80
         MPTST = LXYTST/2                                                 40.80
      ELSE                                                                40.80
         MPTST = LXYTST                                                   40.80
      ENDIF                                                               40.80
      CALL INKEYW ('STA','IJ')                                            40.00
      IF (MCGRD.GT.1 .OR. nverts.GT.0) THEN                               40.80 40.04
        IF (KEYWIS('XY')) THEN
          LOCGRI = .TRUE.
        ELSE IF (KEYWIS('IJ')) THEN
          LOCGRI = .FALSE.
        ELSE
          CALL WRNKEY
        ENDIF
      ELSE                                                                40.04
        CALL MSGERR(3,
     &     'command READ BOT or READ UNSTRUC must precede command TEST')  40.80
      ENDIF                                                               40.04
!
  10  IF (LOCGRI) THEN
        CALL READXY ('X','Y',XP,YP, 'REP', -1.E10, -1.E10)                40.03
        IF (XP.LT.-.9E10) THEN
          LXDMP = 0
          GOTO 60
        ELSE
          IF (OPTG.NE.5) THEN                                             40.80
            CALL CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID, YCGRID, KGRBND)  40.00
            IF (XC.LT.0.) THEN                                            40.31
              IF (XP.GE.XCGMIN .AND. XP.LE.XCGMAX .AND.                   40.31
     &            YP.GE.YCGMIN .AND. YP.LE.YCGMAX ) THEN                  40.31
                 GOTO 50                                                  40.31
              ELSE                                                        40.31
                 GOTO 40                                                  40.31
              END IF                                                      40.31
            END IF                                                        40.31
            LXDMP = NINT(XC) + MXF -1                                     40.31
            LYDMP = NINT(YC) + MYF -1                                     40.31
            IF (ITEST.GE.30) WRITE (PRTEST, 14) XP+XOFFS, YP+YOFFS,       40.03
     &      LXDMP, LYDMP
  14        FORMAT (' test point ', 2F12.2, ' to grid point ', 2I4)
          ELSE                                                            40.80
            CALL SwanFindPoint ( XP, YP, K )                              40.80
            IF ( K.LT.0 ) THEN                                            40.80
              IF (XP.GE.XCGMIN .AND. XP.LE.XCGMAX .AND.                   40.80
     &            YP.GE.YCGMIN .AND. YP.LE.YCGMAX ) THEN                  40.80
                 GOTO 50                                                  40.80
              ELSE                                                        40.80
                 GOTO 40                                                  40.80
              END IF                                                      40.80
            END IF                                                        40.80
            LXDMP = K                                                     40.80
            IF (ITEST.GE.30) WRITE (PRTEST,15) XP+XOFFS,YP+YOFFS,LXDMP    40.80
  15        FORMAT (' test point ', 2F12.2, ' to vertex ', I6)            40.80
          ENDIF                                                           40.80
        ENDIF
      ELSE
        CALL ININTG ('I' , LXDMP, 'REP', -1)                              40.03
        IF (LXDMP .LT. 0) GOTO 60
        IF (OPTG.NE.5) CALL ININTG ('J' , LYDMP, 'REQ',  0)               40.80 40.03
      ENDIF
!
      IF (OPTG.NE.5) THEN                                                 40.80
         IF (LXDMP.GE.0 .AND. LXDMP.LE.MXCGL-1 .AND.                      40.31
     &       LYDMP.GE.0 .AND. LYDMP.LE.MYCGL-1) THEN                      40.31
            LXDMP = LXDMP - MXF + 1                                       40.31
            LYDMP = LYDMP - MYF + 1                                       40.31
            IF (LXDMP.GE.0 .AND. LXDMP.LE.MXC-1 .AND.
     &          LYDMP.GE.0 .AND. LYDMP.LE.MYC-1) THEN
               IF (KGRPNT(LXDMP+1,LYDMP+1) .GT. 1) THEN
                  NPTST = NPTST + 1
                  XYTST(2*NPTST-1) = LXDMP+1
                  XYTST(2*NPTST)   = LYDMP+1
                  GOTO 50
               ENDIF
            ELSE                                                          40.31
               GOTO 50                                                    40.31
            ENDIF                                                         40.31
         ENDIF
      ELSE                                                                40.80
         IF (LXDMP.GE.1 .AND. LXDMP.LE.nverts) THEN                       40.80
            NPTST = NPTST + 1                                             40.80
            XYTST(NPTST) = LXDMP                                          40.80
            GOTO 50                                                       40.80
         ENDIF                                                            40.80
      ENDIF                                                               40.80
!
  40  CALL MSGERR (1, 'test point is not active')                         40.80
      WRITE (PRINTF, *) XP+XOFFS, YP+YOFFS                                40.03
  50  IF (NPTST.LE.MPTST) GOTO 10
      CALL MSGERR (2, 'Too many test points')
!
!     generate output point set 'TESTPNTS'
!
  60  ALLOCATE(OPSTMP)                                                    40.31
      OPSTMP%PSNAME = 'TESTPNTS'                                          40.31
      OPSTMP%PSTYPE = 'P'                                                 40.31
      OPSTMP%MIP    = NPTST                                               40.31
      ALLOCATE(OPSTMP%XP(NPTST))                                          40.31
      ALLOCATE(OPSTMP%YP(NPTST))                                          40.31
      IF (OPTG.NE.5) THEN                                                 40.80
         DO IPTST = 1, NPTST                                              40.80
            LXDMP = XYTST(2*IPTST-1)                                      40.00
            LYDMP = XYTST(2*IPTST)                                        40.00
            OPSTMP%XP(IPTST) = XCGRID(LXDMP,LYDMP)                        40.31
            OPSTMP%YP(IPTST) = YCGRID(LXDMP,LYDMP)                        40.31
         ENDDO                                                            40.80
      ELSE                                                                40.80
         DO IPTST = 1, NPTST                                              40.80
            K = XYTST(IPTST)                                              40.80
            OPSTMP%XP(IPTST) = xcugrd(K)                                  40.80
            OPSTMP%YP(IPTST) = ycugrd(K)                                  40.80
         ENDDO                                                            40.80
      ENDIF                                                               40.80
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
!     open output file for test output of wave parameters                 40.00
!
      CALL INKEYW ('STA', ' ')                                            40.00
      IF (KEYWIS('PAR')) THEN                                             40.00
        CALL INCSTR ('FNAME', FILENM, 'STA', 'SWSRCPA')                   40.00
!       --- append node number to FILENM in case of                       40.30
!           parallel computing                                            40.30
        IF ( PARLL ) THEN                                                 40.30
           ILPOS = INDEX ( FILENM, ' ' )-1                                40.30
           WRITE(FILENM(ILPOS+1:ILPOS+4),99) INODE                        40.30
  99       FORMAT('-',I3.3)                                               40.30
        END IF                                                            40.30
!PUN        FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                     40.95
        IERR = 0                                                          40.00
        CALL FOR (IFPAR, FILENM, 'UF', IERR)                              40.00
        IF (STPNOW()) RETURN                                              34.01
        WRITE (IFPAR, 101) 1                                              40.00
 101    FORMAT ('SWAN', I4, T41,
     &    'Swan standard spectral file, version')                         40.00
        WRITE (IFPAR, 111) VERTXT                                         40.03
 111    FORMAT ('$   Data produced by SWAN version ', A)                  40.03
        WRITE (IFPAR, 113) PROJID, PROJNR                                 40.03
 113    FORMAT ('$   Project: ', A, ';  run number: ', A)
        IF (NSTATM.EQ.1) THEN
          WRITE (IFPAR, 102) 'TIME', 'time-dependent data'
 102      FORMAT (A, T41, A)                                              40.00
          WRITE (IFPAR, 103) ITMOPT, 'time coding option'                 40.03
 103      FORMAT (I6, T41, A)                                             40.00
        ELSE
          WRITE (IFPAR, 102) 'ITER', 'iteration-dependent data'
          WRITE (IFPAR, 103) 0
        ENDIF
        IF (KSPHER.EQ.0) THEN
          WRITE (IFPAR, 102) 'LOCATIONS', 'locations in x-y-space'
        ELSE
          WRITE (IFPAR, 102) 'LONLAT',
     &                       'locations in longitude, latitude'
        ENDIF
        WRITE (IFPAR, 103) NPTST, 'number of locations'
        IF (OPTG.NE.5) THEN                                               40.80
           DO 110 IPTST = 1, NPTST
              LXDMP = XYTST(2*IPTST-1)                                    40.00
              LYDMP = XYTST(2*IPTST)                                      40.00
              WRITE (IFPAR, 106) XCGRID(LXDMP,LYDMP)+XOFFS,
     &                           YCGRID(LXDMP,LYDMP)+YOFFS                40.00
 110       CONTINUE
        ELSE                                                              40.80
           DO IPTST = 1, NPTST                                            40.80
              K = XYTST(IPTST)                                            40.80
              WRITE (IFPAR, 106) xcugrd(K)+XOFFS, ycugrd(K)+YOFFS         40.80
           ENDDO                                                          40.80
        ENDIF                                                             40.80
 106    FORMAT (2(1X,F12.2))
        WRITE (IFPAR, 132) 9                                              40.55 40.00
 132    FORMAT ('QUANT', /, I6, T41, 'number of quantities in table')     40.00
        WRITE (IFPAR, 102) OVSNAM(10), OVLNAM(10)                         40.00
        WRITE (IFPAR, 102) OVUNIT(10), 'unit'                             40.00
        WRITE (IFPAR, 104) OVEXCV(10), 'exception value'                  40.00
 104    FORMAT (F14.6, T41, A)                                            40.00
        WRITE (IFPAR, 102) OVSNAM(28), OVLNAM(28)                         40.41 40.00
        WRITE (IFPAR, 102) OVUNIT(28), 'unit'                             40.41 40.00
        WRITE (IFPAR, 104) OVEXCV(28), 'exception value'                  40.41 40.00
        WRITE (IFPAR, 102) 'Swind',  'wind source term (of var. dens.)'   40.00
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.00
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFPAR, 102) 'Swcap',  'whitecapping dissipation'           40.00
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.00
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFPAR, 102) 'Sfric',  'bottom friction dissipation'        40.00
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.00
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFPAR, 102) 'Svege',  'vegetation dissipation'             40.55
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.55
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.55
        WRITE (IFPAR, 102) 'Ssurf',  'surf breaking dissipation'          40.00
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.00
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFPAR, 102) 'Snl3',   'total absolute 3-wave interaction'  40.13
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.13
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFPAR, 102) 'Snl4',   'total absolute 4-wave interaction'  40.13
        WRITE (IFPAR, 102) 'm2/s',   'unit'                               40.13
        WRITE (IFPAR, 104) OVEXCV(7),'exception value'                    40.00
      ENDIF
!
!     open output file for source terms if requested
!
      CALL INKEYW ('STA', ' ')                                            40.00
      IF (KEYWIS('S1D')) THEN                                             40.00
        CALL INCSTR ('FNAME', FILENM, 'STA', 'SWSRC1D')                   40.00
!       --- append node number to FILENM in case of                       40.30
!           parallel computing                                            40.30
        IF ( PARLL ) THEN                                                 40.30
           ILPOS = INDEX ( FILENM, ' ' )-1                                40.30
           WRITE(FILENM(ILPOS+1:ILPOS+4),99) INODE                        40.30
        END IF                                                            40.30
!PUN        FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                     40.95
        IERR = 0                                                          40.00
        CALL FOR (IFS1D, FILENM, 'UF', IERR)                              40.00
        IF (STPNOW()) RETURN                                              34.01
        WRITE (IFS1D, 101) 1
        WRITE (IFS1D, 111) VERTXT                                         40.03
        WRITE (IFS1D, 113) PROJID, PROJNR                                 40.03
        IF (NSTATM.EQ.1) THEN
          WRITE (IFS1D, 102) 'TIME', 'time-dependent data'
          WRITE (IFS1D, 103) ITMOPT, 'time coding option'                 40.03
        ELSE
          WRITE (IFS1D, 102) 'ITER', 'iteration-dependent data'
          WRITE (IFS1D, 103) 0
        ENDIF
        IF (KSPHER.EQ.0) THEN
          WRITE (IFS1D, 102) 'LOCATIONS', 'locations in x-y-space'
        ELSE
          WRITE (IFS1D, 102) 'LONLAT',
     &                       'locations in longitude, latitude'
        ENDIF
        WRITE (IFS1D, 103) NPTST, 'number of locations'
        IF (OPTG.NE.5) THEN                                               40.80
           DO 210 IPTST = 1, NPTST                                        40.00
              LXDMP = XYTST(2*IPTST-1)                                    40.00
              LYDMP = XYTST(2*IPTST)                                      40.00
              WRITE (IFS1D, 106) XCGRID(LXDMP,LYDMP)+XOFFS,               40.00
     &                           YCGRID(LXDMP,LYDMP)+YOFFS                40.00
 210       CONTINUE                                                       40.00
        ELSE                                                              40.80
           DO IPTST = 1, NPTST                                            40.80
              K = XYTST(IPTST)                                            40.80
              WRITE (IFS1D, 106) xcugrd(K)+XOFFS, ycugrd(K)+YOFFS         40.80
           ENDDO                                                          40.80
        ENDIF                                                             40.80
        IF (ICUR.GT.0) THEN
          WRITE (IFS1D, 102) 'RFREQ', 'relative frequencies in Hz'        40.00
        ELSE
          WRITE (IFS1D, 102) 'AFREQ', 'absolute frequencies in Hz'        40.00
        ENDIF
        WRITE (IFS1D, 103) MSC, 'number of frequencies'                   40.00
        DO 220 IS = 1, MSC                                                40.00
          WRITE (IFS1D, 214) SPCSIG(IS)/PI2                               40.00
 214      FORMAT (F10.4)                                                  40.00
 220    CONTINUE                                                          40.00
        WRITE (IFS1D, 132) 8                                              40.55
        WRITE (IFS1D, 102) 'VaDens', 'variance densities'                 40.00
        WRITE (IFS1D, 102) 'm2/Hz',  'unit'                               40.00
        WRITE (IFS1D, 104) 0.,       'exception value'                    40.00
        WRITE (IFS1D, 102) 'Swind',  'wind source term'                   40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS1D, 102) 'Swcap',  'whitecapping dissipation'           40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS1D, 102) 'Sfric',  'bottom friction dissipation'        40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS1D, 102) 'Svege',  'vegetation dissipation'             40.55
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.55
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.55
        WRITE (IFS1D, 102) 'Ssurf',  'surf breaking dissipation'          40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS1D, 102) 'Snl3',   'triad interactions'                 40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS1D, 102) 'Snl4',   'quadruplet interactions'            40.00
        WRITE (IFS1D, 102) 'm2',     'unit'                               40.00
        WRITE (IFS1D, 104) OVEXCV(7),'exception value'                    40.00
      ENDIF
!     2D source terms
      CALL INKEYW ('STA', ' ')                                            40.00
      IF (KEYWIS('S2D')) THEN                                             40.00
        CALL INCSTR ('FNAME', FILENM, 'STA', 'SWSRC2D')                   40.00
!       --- append node number to FILENM in case of                       40.30
!           parallel computing                                            40.30
        IF ( PARLL ) THEN                                                 40.30
           ILPOS = INDEX ( FILENM, ' ' )-1                                40.30
           WRITE(FILENM(ILPOS+1:ILPOS+4),99) INODE                        40.30
        END IF                                                            40.30
!PUN        FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)                     40.95
        IERR = 0                                                          40.00
        CALL FOR (IFS2D, FILENM, 'UF', IERR)                              40.00
        IF (STPNOW()) RETURN                                              34.01
        WRITE (IFS2D, 101) 1
        WRITE (IFS2D, 111) VERTXT                                         40.03
        WRITE (IFS2D, 113) PROJID, PROJNR                                 40.03
        IF (NSTATM.EQ.1) THEN
          WRITE (IFS2D, 102) 'TIME', 'time-dependent data'
          WRITE (IFS2D, 103) ITMOPT, 'time coding option'                 40.03
        ELSE
          WRITE (IFS2D, 102) 'ITER', 'iteration-dependent data'
          WRITE (IFS2D, 103) 0
        ENDIF
        IF (KSPHER.EQ.0) THEN
          WRITE (IFS2D, 102) 'LOCATIONS', 'locations in x-y-space'
        ELSE
          WRITE (IFS2D, 102) 'LONLAT',
     &                       'locations in longitude, latitude'
        ENDIF
        WRITE (IFS2D, 103) NPTST, 'number of locations'
        IF (OPTG.NE.5) THEN                                               40.80
           DO 310 IPTST = 1, NPTST                                        40.00
              LXDMP = XYTST(2*IPTST-1)                                    40.00
              LYDMP = XYTST(2*IPTST)                                      40.00
              WRITE (IFS2D, 106) XCGRID(LXDMP,LYDMP)+XOFFS,               40.00
     &                           YCGRID(LXDMP,LYDMP)+YOFFS                40.00
 310       CONTINUE                                                       40.00
        ELSE                                                              40.80
           DO IPTST = 1, NPTST                                            40.80
              K = XYTST(IPTST)                                            40.80
              WRITE (IFS2D, 106) xcugrd(K)+XOFFS, ycugrd(K)+YOFFS         40.80
           ENDDO                                                          40.80
        ENDIF                                                             40.80
        IF (ICUR.GT.0) THEN
          WRITE (IFS2D, 102) 'RFREQ', 'relative frequencies in Hz'        40.00
        ELSE
          WRITE (IFS2D, 102) 'AFREQ', 'absolute frequencies in Hz'        40.00
        ENDIF
        WRITE (IFS2D, 103) MSC, 'number of frequencies'                   40.00
        DO 320 IS = 1, MSC                                                40.00
          WRITE (IFS2D, 214) SPCSIG(IS)/PI2                               40.00
 320    CONTINUE                                                          40.00
!       full 2-D spectrum
        WRITE (IFS2D, 102) 'CDIR',
     &                      'spectral Cartesian directions in degr'       40.00
        WRITE (IFS2D, 103) MDC, 'number of directions'
        DO 330 ID = 1, MDC                                                40.00
          WRITE (IFS2D, 324) SPCDIR(ID,1)*180./PI                         30.82
 324      FORMAT (F10.4)                                                  40.00
 330    CONTINUE                                                          40.00
        WRITE (IFS2D, 132) 8                                              40.55
        WRITE (IFS2D, 102) 'VaDens', 'variance densities'                 40.00
        WRITE (IFS2D, 102) 'm2/Hz/degr', 'unit'                           40.00
        WRITE (IFS2D, 104) 0.,       'exception value'                    40.00
        WRITE (IFS2D, 102) 'Swind',  'wind source term'                   40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS2D, 102) 'Swcap',  'whitecapping dissipation'           40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS2D, 102) 'Sfric',  'bottom friction dissipation'        40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS2D, 102) 'Svege',  'vegetation dissipation'             40.55
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.55
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.55
        WRITE (IFS2D, 102) 'Ssurf',  'surf breaking dissipation'          40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS2D, 102) 'Snl3',   'triad interactions'                 40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
        WRITE (IFS2D, 102) 'Snl4',   'quadruplet interactions'            40.00
        WRITE (IFS2D, 102) 'm2/degr','unit'                               40.00
        WRITE (IFS2D, 104) OVEXCV(7),'exception value'                    40.00
      ENDIF
      RETURN
!     end of subroutine RETSTP
      END
