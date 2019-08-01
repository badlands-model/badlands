!     Last change:  NB   04 Jan 2001   12:00 pm
!
!          OCEAN PACK - Installation dependent subroutines
!
!*****************************************************************
!                                                                *
      SUBROUTINE OCPINI (INIFIL, INPUTFILE, INFONAME, LREAD, INERR)
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
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
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     34.01: IJsbrand Haagsma
!     40.00, 40.03: Nico Booij
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.95: Marcel Zijlema
!
!  1. Updates
!
!     10.02, July 94: New argument INIFIL
!                     Check on validity period now uses OCDTIM
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.82, Nov. 98: Introduced recordlength of 1000 for file PRINT to
!                     avoid error-messages on the Cray-J90 and SGI Origin 200
!     34.01, Feb. 99: Changed STOP statements for MSGERR(4,'message')
!                     calls
!     34.01, Feb. 99: Opens a file 'screen' when unitnr in swaninit<>6
!     40.00, Feb. 99: Directory separation characters included in init file
!                     these characters are used in subr FOR
!     40.03, May  00: backslash replaced by CHAR(92) because of problems on Linux
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing HPGL-functionality
!     40.41, Sep. 04: includes speed processors in initialisation file
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.95, Jun. 08: parallelization of unSWAN using MESSENGER of ADCIRC
!
!  2. Purpose
!
!     subroutine initialises a number of common variables
!     opens standard input and output files, if necessary
!
!  4. Argument variables
!
!     INERR : output     Number of the initialisation error
!
      INTEGER INERR
!
!     INIFIL  inp  char  name of initialisation file
!     LREAD   inp  log   if True: command input file must be opened
!                        and command reading must be initialised
!
      LOGICAL    LREAD, FILEXI
      CHARACTER  INPFIL *40, OUTFIL *40, INIFIL *(*), TSTFIL *40,
     &           PLTOPT *4, TIMSTR *24, INPUTFILE*40, INFONAME*40,
     &           INPFO *40, OUTFO *40, TSTFO *40, TXT*120
      INTEGER    PRCTIM(6), INIVER, INIVEF
      INTEGER    PFROPT                                                   40.31
      REAL       PLPARM(10)                                               40.31
      INTEGER    NUMM(10)                                                 40.41
      LOGICAL    STPNOW
      DATA PRCTIM /0,0,0,0,0,0/
!
!     version of initialisation file                                      10.16
      INIVER = 4                                                          40.41 40.31
      INIVEF = -1                                                         40.00
      INERR  = 0                                                          34.01
!
      IF (PARLL) THEN
         IF (.NOT.ALLOCATED(IWEIG)) ALLOCATE(IWEIG(NPROC))
         IWEIG = 100
      END IF
!
!     see whether initialisation file exists
!
      INQUIRE (FILE=INIFIL, EXIST=FILEXI)
      IF (FILEXI) THEN
!
!       read initialisation file
!
        OPEN (11, FILE=INIFIL, STATUS='OLD',                              40.41
!CVIS     &        SHARED,                                                     40.41
     &        ERR=950)                                                    40.41
        READ (11, *, ERR=930, END=930) INIVEF                             10.16
        IF (INIVEF.GT.INIVER .OR. INIVEF.LE.0) GOTO 935                   40.00
        READ (11, 120, ERR=930, END=930) INST
        READ (11, *,   ERR=930, END=930) INPUTF
        READ (11, 120, ERR=930, END=930) INPFIL
        READ (11, *,   ERR=930, END=930) PRINTF
        READ (11, 120, ERR=930, END=930) OUTFIL
        READ (11, *,   ERR=930, END=930) PRTEST
        READ (11, 120, ERR=930, END=930) TSTFIL
        READ (11, *,   ERR=930, END=930) SCREEN
        READ (11, *,   ERR=930, END=930) IUNMAX
        READ (11, 130, ERR=930, END=930) COMID
        READ (11, 130, ERR=930, END=930) TABC
        IF (INIVEF.GE.2) THEN                                             40.00
          READ (11, 130, ERR=930, END=930) DIRCH1
          READ (11, 130, ERR=930, END=930) DIRCH2
        ELSE
!DOS          DIRCH1 = CHAR(47)                                               40.41 40.03
!DOS          DIRCH2 = CHAR(92)                                               40.03
          DIRCH1 = CHAR(92)                                               40.41 40.30
          DIRCH2 = CHAR(47)                                               40.30
        ENDIF
        IF (INIVEF.LT.3) THEN
           READ (11, 140, ERR=930, END=930) PLTOPT
           READ (11, *,   ERR=930, END=930) NPLP
           READ (11, *,   ERR=930, END=930) (PLPARM(II),II=1,NPLP)
           READ (11, *,   ERR=930, END=930) PFROPT
        END IF
        READ (11, *,   ERR=930, END=930) ITMOPT
        IF (INIVEF.GT.3) THEN
           IF (PARLL) THEN
              DO JJ = 1, NPROC
                 READ (11, 145, ERR=150, END=150) IW, TXT
                 CALL TXPBLA(TXT,IF,IL)
                 IPOS = 0
                 DO II = IF, IL
                    K = ICHAR(TXT(II:II))
                    IF ( K.GE.48 .AND. K.LE.57 ) THEN
                       IPOS = IPOS + 1
                       NUMM(IPOS) = K - 48
                    END IF
                 END DO
                 J = 0
                 DO II = 1, IPOS
                    J = J + NUMM(II)*10**(IPOS-II)
                 END DO
                 IWEIG(J) = IW
              END DO
           END IF
        END IF
 120    FORMAT (A40)
 130    FORMAT (A1)
 140    FORMAT (A4)
 145    FORMAT (I5,A)
 150    CLOSE (11)
      ELSE
!
!       REFERENCE NUMBERS AND NAMES OF STANDARD FILES
!
        INPUTF = 3
!        INPFIL = 'INPUT'
        INPFIL = INPUTFILE
        PRINTF = 4
        OUTFIL = INFONAME
!       unit ref. numbers for output to screen and to separate
!       test print file:
        PRTEST = PRINTF
        TSTFIL = '    '
        SCREEN = 6
        IUNMAX = 99
!/SGI        IUNMAX = 199
!       TABC is the Tab character (interpreted as blank in command reading)
        TABC = CHAR(9)                                                    40.69
!       COMID is the comment identifier (usually $)
        COMID  = '$'
!       DIRCH1 is directory separation character as appears in input file
!       DIRCH2 is directory separation character replacing DIRCH1         40.00
!DOS        DIRCH1 =  CHAR(47)                                                40.41 40.03
!DOS        DIRCH2 =  CHAR(92)                                                40.03
        DIRCH1 =  CHAR(92)                                                40.41 40.30
        DIRCH2 =  CHAR(47)                                                40.30
!       INST = name of institute, max. 40 characters
        INST = 'Delft University of Technology'
        ITMOPT = 1
!
      ENDIF
      INPFO = INPFIL                                                      40.95
      OUTFO = OUTFIL                                                      40.30
      TSTFO = TSTFIL                                                      40.30

!     --- append node number to OUTFIL and TSTFIL                         40.30
!         in case of parallel computing                                   40.30

      IF (PARLL) THEN                                                     40.30
         ILPOS = INDEX ( OUTFIL, ' ' )-1
         IF (ILPOS.GT.0) THEN
            WRITE(OUTFIL(ILPOS+1:ILPOS+4),180) INODE
         ELSE
            GOTO 920
         END IF
         ILPOS = INDEX ( TSTFIL, ' ' )-1
         IF (ILPOS.GT.0) THEN
            WRITE(TSTFIL(ILPOS+1:ILPOS+4),180) INODE
         END IF
 180     FORMAT('-',I3.3)
         CALL SWSYNC
         IF (STPNOW()) RETURN
      END IF

!PUN      IAMMASTER = MYPROC.EQ.0                                             40.95
!PUN      IF ( MNPROC==1 ) THEN                                               40.95
!PUN         INPUTDIR = '.'                                                   40.95
!PUN         LOCALDIR = '.'                                                   40.95
!PUN      ELSE                                                                40.95
!PUN         CALL MAKE_DIRNAME()                                              40.95
!PUN      ENDIF                                                               40.95
!PUN      INPFIL = TRIM(INPUTDIR)//DIRCH2//TRIM(INPFIL)                       40.95
!PUN      OUTFIL = TRIM(LOCALDIR)//DIRCH2//TRIM(OUTFIL)                       40.95
!PUN      TSTFIL = TRIM(LOCALDIR)//DIRCH2//TRIM(TSTFIL)                       40.95

!      IF (INIVEF.LT.INIVER) THEN
!
!       write initialisation file
!
!        OPEN  (12, FILE=INIFIL, STATUS='UNKNOWN', FORM='FORMATTED',       40.00
!     &           ERR=950)
!        WRITE (12, 210) INIVER, 'version of initialisation file'          10.16
!        WRITE (12, 220) INST,   'name of institute'
!        WRITE (12, 210) INPUTF, 'command file ref. number'
!        WRITE (12, 220) INPFO,  'command file name'
!        WRITE (12, 210) PRINTF, 'print file ref. number'
!        WRITE (12, 220) OUTFO,  'print file name'                         40.30
!        WRITE (12, 210) PRTEST, 'test file ref. number'
!        WRITE (12, 220) TSTFO,  'test file name'                          40.30
!        WRITE (12, 210) SCREEN, 'screen ref. number'
!        WRITE (12, 210) IUNMAX, 'highest file ref. number'
!        WRITE (12, 230) COMID,  'comment identifier'
!        WRITE (12, 230) TABC,   'TAB character'
!        WRITE (12, 230) DIRCH1, 'dir sep char in input file'
!        WRITE (12, 230) DIRCH2, 'dir sep char replacing previous one'     40.00
!        WRITE (12, 210) ITMOPT, 'default time coding option'
!        IF (PARLL) THEN
!           DO II = 1, NPROC
!              WRITE (12, 240) IWEIG(II), 'speed of processor ',II
!           END DO
!        END IF
!        CLOSE (12)
!  210   FORMAT (I5, T41, A)
!  220   FORMAT (A40, A)
!  230   FORMAT (A1, T41, A)
!  240   FORMAT (I5, T41, A19, I3)
!      ENDIF
!
      IUNMIN = 0
      FUNLO = 21
!/SGI      FUNLO = 103
      FUNHI = IUNMAX
!
      CALL OCDTIM (PRCTIM)
!
!     initialise command reader
!
      IF (OUTFIL.NE.'    ') THEN
!       WRITE (*,*) ' Open print file ', PRINTF, OUTFIL
        OPEN (UNIT=PRINTF, FILE=OUTFIL, STATUS='UNKNOWN',
     &    FORM='FORMATTED',                                               30.82
!/Cray     &    RECL=2000,                                                      30.82
!/SGI     &    RECL=2000,                                                      30.82
     &    ERR=920)                                                        30.82
!       WRITE (*,*) ' Print file opened ', PRINTF, OUTFIL
        CALL DTTIST (ITMOPT, TIMSTR, PRCTIM)
        WRITE (PRINTF, 12) TIMSTR
   12   FORMAT ('1',//,20X, 'Execution started at ',A, //)
      ENDIF
      IF (PRTEST.NE.PRINTF) OPEN (UNIT=PRTEST, FILE=TSTFIL, ERR=922)
      IF (SCREEN.NE.6) OPEN(UNIT=SCREEN, FILE='screen',ERR=960)           34.01
      IF (LREAD) THEN
        IF (INPFIL.NE.'    ')
     &  OPEN (UNIT=INPUTF, FILE=INPFIL, STATUS='OLD',                     40.41
!CVIS     &        SHARED,                                                     40.41
     &        ERR=910)                                                    40.41
        CALL RDINIT
      ENDIF
!
      RETURN
!
  910 CALL MSGERR(4,'Input file missing')                                 34.01
      RETURN                                                              34.01
!
  920 INERR=920                                                           34.01
      IF ( IAMMASTER )                                                    40.95 40.31
     &   WRITE(*,*) 'Cannot open PRINT file '                             40.31
      RETURN                                                              34.01
!
  922 CALL MSGERR(4,'Cannot open test file: '//TSTFIL)                    34.01
      RETURN
!
  930 INERR=930                                                           34.01
      IF ( IAMMASTER )                                                    40.95 40.31
     &   WRITE(*,*) 'Error reading initialisation file '                  40.31
      RETURN                                                              34.01
!
  935 INERR=935                                                           34.01
      IF ( IAMMASTER )                                                    40.95 40.31
     &   WRITE(*,*) 'Incorrect version of initialisation file '           40.31
      RETURN                                                              34.01
!
  950 INERR=950                                                           34.01
      IF ( IAMMASTER )                                                    40.95 40.31
     &   WRITE(*,*) 'Error opening initialisation file '                  40.31
      RETURN                                                              34.01
!
  960 CALL MSGERR(4,'Error opening output file: screen')                  34.01
      RETURN                                                              34.01
!
      END
!*****************************************************************
!                                                                *
      SUBROUTINE OCDTIM (PRCTIM)
!                                                                *
!*****************************************************************
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
!     30.07
!     30.70: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!
!  1. Updates
!
!     30.07, Oct. 95: option DEC added
!     30.70, Sep. 97: adaptation in view of year 2000
!     30.82, Mar. 99: Adapted to Fortran 90 standard
!     40.02, Sep. 00: Removed all platform dependent Fortran 77 statements
!
!  2. PURPOSE
!
!       get time of processing, using processor dependent routines
!
!  3. PARAMETER LIST
!
!       PRCTIM  outp   int   time array: elements: year, month, day,
!                            hour, minute, second
!
!  4. SUBROUTINES USED
!
!       GETDAT, GETTIM or other
!
!  5. ERROR MESSAGES
!
!       ----
!
!  6. REMARKS
!
!       This function uses a processor dependent subroutines (GETDAT,
!       GETTIM), therefore adaptations are necessary when compiled
!       at a different computer system environment.
!
!  7. STRUCTURE
!
!       ---------------------------------------------------------
!       Call DATE and TIME routines (system dependent)
!       decode YEAR, MONTH and DAY
!       assemble DATE string
!       decode HOUR, MINUTE, SECOND
!       assemble TIME string
!       ---------------------------------------------------------
!
!  8. SOURCE TEXT
!
      INTEGER PRCTIM(6)
!
!     Call DATE and TIME routines
!
!     --------Fortran 90 date-time routines --------
!
      CHARACTER TIMSTR *24, CDUMMY *5                                     30.82
      INTEGER   IDUMMY(8)                                                 30.82
!
      CALL DATE_AND_TIME (TIMSTR(1:8), TIMSTR(10:20), CDUMMY, IDUMMY)     30.70
      CALL DTSTTI (1, TIMSTR, PRCTIM)                                     30.70

      RETURN

      END SUBROUTINE OCDTIM
!*****************************************************************
!                                                                *
      SUBROUTINE DTSTTI (IOPT, TIMSTR, DTTIME)
!                                                                *
!*****************************************************************
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
!     Updates
!
!       ver 30.70, Sep 1997 by N.Booij: adaptation in view year 2000
!       40.89, Nico Booij: integer Jumpyear introduced to prevent porblems
!                          with 2-digit year code
!
!     Function:
!
!       transform time string into integer time array
!
!     Argument list:
!
!       IOPT    input  int   option number
!                            1: ISO notation   19870530.153000
!                            2: (HP compiler): 30-May-87 15:30:00
!                            3: (old Lahey)    05/30/87 15:30:00          30.70
!                            4:                         15:30:00
!                            5:                87/05/30 15:30:00
!                            6: WAM            8705301530
!
!       TIMSTR  input  char  time string
!       DTTIME  outp   int   time array: elements: year, month, day,
!                            hour, minute, second
!
!    Remarks:
!     Options can be added by the user
!     existing options should not be changed
!
!     Source:
!
      INTEGER :: JUMPYEAR = 30     ! Used in interpretation of 2-digit year code
      INTEGER    IOPT, DTTIME(6)
      CHARACTER  TIMSTR *24, MONC(12) *3, MONCI *3
      DATA MONC /'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
     &           'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/
!
      IF (IOPT.EQ.1) THEN
        READ (TIMSTR, '(I4,I2,I2,1X,3I2)', ERR=98) (DTTIME(II), II=1,6)
      ELSE IF (IOPT.EQ.2) THEN
        READ (TIMSTR, '(I2,1X,A3,1X,I2,3(1X,I2))', ERR=98)
     &  DTTIME(3), MONCI, DTTIME(1), (DTTIME(II), II=4,6)
        IF (DTTIME(1).LT.JUMPYEAR) THEN                                   40.89 30.70
          DTTIME(1) = 2000 + DTTIME(1)                                    30.70
        ELSE
          DTTIME(1) = 1900 + DTTIME(1)
        ENDIF
        DTTIME(2) = 0
        DO 20 IMM = 1, 12                                                  !!
          CALL UPCASE (MONCI)
          IF (MONCI.NE.MONC(IMM)) GOTO 20
          DTTIME(2) = IMM
          GOTO 90
  20    CONTINUE                                                           !!
        CALL MSGERR (2, 'incorrect month string: '//MONCI)                40.00
      ELSE IF (IOPT.EQ.3) THEN
        READ (TIMSTR, '(I2,5(1X,I2))', ERR=98)
     &  DTTIME(2), DTTIME(3), DTTIME(1), (DTTIME(II), II=4,6)
        IF (DTTIME(1).LT.JUMPYEAR) THEN                                   40.89 30.70
          DTTIME(1) = 2000 + DTTIME(1)                                    30.70
        ELSE
          DTTIME(1) = 1900 + DTTIME(1)
        ENDIF
      ELSE IF (IOPT.EQ.4) THEN
        READ (TIMSTR, '(I2,2(1X,I2))', ERR=98) (DTTIME(II), II=4,6)
        DO 40 II = 1, 3
          DTTIME(II) = 0
  40    CONTINUE
      ELSE IF (IOPT.EQ.5) THEN                                            30.11
        READ (TIMSTR, '(I2,5(1X,I2))', ERR=98) (DTTIME(II), II=1,6)       30.11
        IF (DTTIME(1).LT.JUMPYEAR) THEN                                   40.89 30.70
          DTTIME(1) = 2000 + DTTIME(1)                                    30.70
        ELSE
          DTTIME(1) = 1900 + DTTIME(1)
        ENDIF
      ELSE IF (IOPT.EQ.6) THEN                                            30.11
        READ (TIMSTR, '(5I2)', ERR=98) (DTTIME(II), II=1,5)
        DTTIME(6) = 0.
        IF (DTTIME(1).LT.JUMPYEAR) THEN                                   40.89 30.70
          DTTIME(1) = 2000 + DTTIME(1)                                    30.70
        ELSE
          DTTIME(1) = 1900 + DTTIME(1)
        ENDIF
      ELSE
        CALL MSGERR (2, 'wrong time coding option in subroutine DTSTTI')
      ENDIF
  90  RETURN
  98  CALL MSGERR (2, 'time string unreadable: '//TIMSTR)
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE DTTIST (IOPT, TIMSTR, DTTIME)
!                                                                *
!*****************************************************************
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
!     Updates
!
!       ver 30.70, Sep 1997 by N.Booij: adaptation in view year 2000
!
!     Function:
!
!       transform integer time array into time string
!
!     Argument list:
!
!       IOPT    input  int   option number (see subr. DTSTTI)
!       TIMSTR  outp   char  time string
!       DTTIME  input  int   time array: elements: year, month, day,
!                            hour, minute, second
!
!     Source:
!
      INTEGER    IOPT, DTTIME(6)
      CHARACTER  TIMSTR *24, MONC(12) *3
      DATA MONC /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
     &           'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
      TIMSTR = '    '
      IF (IOPT.EQ.1) THEN
        WRITE (TIMSTR, 12) (DTTIME(II), II=1,6)
  12    FORMAT (I4,I2,I2,'.',3I2)
        LTS = 15
      ELSE IF (IOPT.EQ.2) THEN
        IF (DTTIME(1).GE.2000) THEN                                       30.70
          DTTIME(1) = DTTIME(1) - 2000                                    30.70
        ELSE
          DTTIME(1) = DTTIME(1) - 1900                                    30.70
        ENDIF
        WRITE (TIMSTR, 22) DTTIME(3), MONC(DTTIME(2)), DTTIME(1),         30.70
     &  (DTTIME(II), II=4,6)
  22    FORMAT (I2,'-',A3,'-',I2,'.',I2,':',I2,':',I2)
        LTS = 18
      ELSE IF (IOPT.EQ.3) THEN
        IF (DTTIME(1).GE.2000) THEN                                       30.70
          DTTIME(1) = DTTIME(1) - 2000                                    30.70
        ELSE
          DTTIME(1) = DTTIME(1) - 1900                                    30.70
        ENDIF
        WRITE (TIMSTR, 32)
     &  DTTIME(2), DTTIME(3), DTTIME(1), (DTTIME(II), II=4,6)             30.70
  32    FORMAT (I2,'/',I2,'/',I2,'.',I2,':',I2,':',I2)
        LTS = 17
      ELSE IF (IOPT.EQ.4) THEN
        WRITE (TIMSTR, 42) (DTTIME(II), II=4,6)
  42    FORMAT (I2,':',I2,':',I2)
        LTS = 8
      ELSE IF (IOPT.EQ.5) THEN                                            30.11
        IF (DTTIME(1).GE.2000) THEN                                       30.70
          DTTIME(1) = DTTIME(1) - 2000                                    30.70
        ELSE
          DTTIME(1) = DTTIME(1) - 1900                                    30.70
        ENDIF
        WRITE (TIMSTR, 52) (DTTIME(II), II= 1,6)                          30.70
  52    FORMAT (I2,'/',I2,'/',I2,'.',I2,':',I2,':',I2)
        LTS = 17
      ELSE IF (IOPT.EQ.6) THEN                                            30.11
        IF (DTTIME(1).GE.2000) THEN                                       30.70
          DTTIME(1) = DTTIME(1) - 2000                                    30.70
        ELSE
          DTTIME(1) = DTTIME(1) - 1900                                    30.70
        ENDIF
        WRITE (TIMSTR, 62) (DTTIME(II), II=1,5)                           30.70
  62    FORMAT (5I2)
        LTS = 10
      ELSE
        CALL MSGERR (2, 'wrong time coding option in subroutine DTTIST')
      ENDIF
!
      DO 94 IC = 1, LTS
        IF (TIMSTR(IC:IC).EQ.' ') TIMSTR(IC:IC) = '0'
  94  CONTINUE
!
      RETURN
      END
