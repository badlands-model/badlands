!     Ocean Pack miscellaneous routines
!
!     real function DTTIME
!     subroutine    DTINTI
!     subroutine    DTRETI
!     char function DTTIWR
!     REPARM
!     INAR2D
!     STRACE
!     MSGERR
!     TABHED
!     FOR
!     logical function EQREAL /* Checks whether REAL1 is appr.            30.72
!                                equal to REAL2                 */        30.72
!     LSPLIT            /* splits an input line into data items */        40.00
!     BUGFIX                                                              40.03
!     COPYCH (copied from file OCPDPN)                                    40.31
!
!*******************************************************************
!                                                                  *
      REAL FUNCTION DTTIME (INTTIM)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!     30.74: IJsbrand Haagsma (Include version)
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!      9705, May  97: month number is checked
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     DTTIME gives time in seconds from a reference day
!            it also initialises the reference day
!
!  3. Method
!
!     every fourth year is a leap-year, but not the century-years, however
!     also leap-years are: year 0, 1000, 2000 etc.
!     1 jan of year 0 is daynumber 1.
!
!  4. Argument variables
!
!     INTTIM(1): year
!           (2): month
!           (3): day
!           (4): hour
!           (5): minute
!           (6): second
!
      INTEGER INTTIM(6)                                                   30.74
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IDYMON : number of days of each month (February counts as 28 days)
!     IYEAR  : number of years after substacking the centuries
!     IYRM1  : ??
!     IDNOW  : ??
!     I      : ??
!     II     : ??
!
      INTEGER IDYMON(12), IYEAR, IYRM1, IDNOW, I, II
!
!     LEAPYR : Whether year in INTTIM(1) is a leapyear
!     LOGREF : ??
!
      LOGICAL LEAPYR, LOGREF
!
!     REFDAY  day number of the reference day; the reference time is 0:00
!            of the reference day; the first day entered is used as
!             reference day.
!
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE LOGREF, IDYMON
      DATA IDYMON /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      DATA LOGREF /.FALSE./
!
      IYEAR = INTTIM(1)
      IYRM1 = IYEAR-1
      LEAPYR=(MOD(IYEAR,4).EQ.0.AND.MOD(IYEAR,100).NE.0).OR.
     &        MOD(IYEAR,1000).EQ.0
      IDNOW=0
      IF (INTTIM(2).GT.12) THEN                                           9705
        WRITE (PRINTF, 8) INTTIM(2), (INTTIM(II), II=1,6)                 9705
   8    FORMAT (' erroneous month ', I2, ' in date/time ', 6I4)           9705
      ELSE IF (INTTIM(2).GT.1) THEN                                       9705
        DO 10 I = 1,INTTIM(2)-1
          IDNOW=IDNOW+IDYMON(I)
  10    CONTINUE
      ENDIF                                                               9705
      IDNOW=IDNOW+INTTIM(3)
      IF (LEAPYR.AND.INTTIM(2).GT.2) IDNOW=IDNOW+1
      IDNOW = IDNOW + IYEAR*365 + IYRM1/4 - IYRM1/100 + IYRM1/1000 + 1
      IF (IYEAR.EQ.0) IDNOW=IDNOW-1
      IF (.NOT.LOGREF) THEN
        REFDAY = IDNOW
        LOGREF = .TRUE.
        DTTIME = 0.
      ELSE
        DTTIME = REAL(IDNOW-REFDAY) * 24.*3600.
      ENDIF
      DTTIME = DTTIME + 3600.*REAL(INTTIM(4)) + 60.*REAL(INTTIM(5)) +
     &                  REAL(INTTIM(6))
      RETURN
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE DTINTI (TIMESC, INTTIM)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        40.41
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
!     30.74: IJsbrand Haagsma (Include version)
!     30.70: Nico Booij (small change)
!
!  1. Updates
!
!      9705, May  97: month number is checked
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.70, Jan. 98: small change in interpretation of time in sec
!
!  2. Purpose
!
!     DTINTI calculates integer time array INTTIM from time in seconds
!            from given reference day REFDAY
!
!  3. Method
!
!     every fourth year is a leap-year, but not the century-years, however
!     also leap-years are: year 0, 1000, 2000 etc.
!     1 jan of year 0 is daynumber 1.
!
!  4. Argument variables
!
!     INTTIM(1): year
!           (2): month
!           (3): day
!           (4): hour
!           (5): minute
!           (6): second
!
      INTEGER INTTIM(6)
!
!     TIMESC : input  time in seconds from given reference day REFDAY
!
      REAL TIMESC
!
!  5. PARAMETER VARIABLES
!
!     IDAYYR : number of days in 'normal' year (no leap-year)
!     IDYCEN : number of days in a century
!     IDYMIL : number of days in a millenium (1000 years)
!     IFOUR  : number of days in 4 year with 1 leap-year
!
      INTEGER IDAYYR, IDYCEN, IDYMIL, IFOUR
!
      PARAMETER (IDAYYR = 365)
      PARAMETER (IDYMIL = IDAYYR*1000+1000/4-1000/100+1)
      PARAMETER (IDYCEN = IDAYYR*100+100/4-1)
      PARAMETER (IFOUR  = 4*IDAYYR+1)
!
!  6. LOCAL VARIABLES
!
!     I4     : number of blocks of four years after subtraction of the
!              millenia and the centuries
!     ICEN   : number of centuries after subtracking the millenia
!     IDYMN  : day of the month
!     IDYMON : number of days of each month (February counts as 28 days)
!     IDYNOW : local daynumber
!     IMIL   : number of millenia in julday-1 days
!     IMN    : month counter
!     IYR    : remaining number of years
!     IYEAR  : number of years after substacking the centuries
!     NDAY   : number of days since reference day
!     NOWDAY : reference day
!
      INTEGER I4, ICEN, IDYMN, IDYMON(12), IDYNOW, IMIL, IMN, IYR, IYEAR
     &, NDAY, NOWDAY
!
!     TT     : time in seconds since begin of the same day
!
      REAL    TT
!
!     LEAPYR : logical for yes or no leap-year
!
      LOGICAL LEAPYR
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE IDYMON
      DATA IDYMON /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
      NDAY = INT((TIMESC+0.4)/(24*3600))                                  30.70
  10  TT   = TIMESC - REAL(NDAY)*24.*3600.
      IF (TT.LT.-0.4) THEN                                                30.70
        NDAY = NDAY - 1
        GOTO 10
      ENDIF
      NOWDAY = REFDAY + NDAY
!
!        get year
!
      IDYNOW = NOWDAY-1
      IMIL   = IDYNOW/IDYMIL
      IDYNOW = IDYNOW-IMIL*IDYMIL
      ICEN   = (IDYNOW-(IDYCEN+1))/IDYCEN+1
      IF (IDYNOW-(IDYCEN+1).LT.0) ICEN=0
      IF (ICEN.EQ.0) THEN
        I4     = IDYNOW/IFOUR
        IDYNOW = IDYNOW-I4*IFOUR
      ELSE
        IDYNOW = IDYNOW-(IDYCEN+1)-(ICEN-1)*IDYCEN
        I4     = (IDYNOW-(IFOUR-1))/IFOUR+1
        IF(IDYNOW-(IFOUR-1).LT.0) I4=0
        IF(I4.GT.0) IDYNOW=IDYNOW-(IFOUR-1)-(I4-1)*IFOUR
      END IF
      IYR   = (IDYNOW-(IDAYYR+1))/IDAYYR+1
      IF(IDYNOW-(IDAYYR+1).LT.0) IYR=0
      IYEAR = 1000*IMIL + 100*ICEN + 4*I4 + IYR
!
!        get month and day
!
      LEAPYR = (MOD(IYEAR,4).EQ.0.AND.MOD(IYEAR,100).NE.0).OR.
     &          MOD(IYEAR,1000).EQ.0
      IF (IYR.GT.0) IDYNOW=IDYNOW-(IDAYYR+1)-(IYR-1)*IDAYYR
      IDYNOW = IDYNOW+1
      DO 30 IMN = 1, 12
        IDYMN=IDYMON(IMN)
        IF(LEAPYR.AND.IMN.EQ.2) IDYMN=IDYMN+1
        IF(IDYNOW.LE.IDYMN) GOTO 40
        IDYNOW=IDYNOW-IDYMN
  30  CONTINUE
  40  INTTIM(2) = IMN
      INTTIM(3) = IDYNOW
      INTTIM(1) = IYEAR
!
!        get time of day
!
      INTTIM(4) = INT(TT/3600.)
      TT        = TT - 3600.*REAL(INTTIM(4))
      INTTIM(5) = INT(TT/60.)
      TT        = TT - 60.*REAL(INTTIM(5))
      INTTIM(6) = INT(TT)
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE DTRETI (TSTRNG, IOPT, TIMESC)
!                                                                *
!*****************************************************************
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
!  0. AUTHORS
!
!  1. UPDATES
!
!  2. PURPOSE
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IOPT   : input    option number
!
      INTEGER IOPT
!
!     TIMESC : output   time in seconds from given reference day REFDAY
!
      REAL    TIMESC
!
!     TSTRNG : input    time string
!
      CHARACTER  TSTRNG *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     ITIME  : ??
!
      INTEGER ITIME(6)
!
!     DTTIME : Gives time in seconds from a reference day it also initialises the
!              reference day
!
      REAL    DTTIME
!
!  8. SUBROUTINE USED
!
!     DTSTTI   (installation dependent subroutines)
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      CALL DTSTTI (IOPT, TSTRNG, ITIME)
      TIMESC = DTTIME (ITIME)
      RETURN
      END
!*****************************************************************
!                                                                *
      CHARACTER *18 FUNCTION DTTIWR (IOPT, TIMESC)                        30.00
!                                                                *
!*****************************************************************
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
!  0. AUTHORS
!
!     40.02: IJsbrand Haagsma
!
!  1. UPDATES
!
!     30.05: New subroutine
!     40.02, Oct. 00: Made length of TSTRNG equal to TIMESTR
!
!  2. PURPOSE
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IOPT   : input    time coding option number
!
      INTEGER    IOPT
!
!     TIMESC : output   time in seconds from given reference day REFDAY
!
      REAL       TIMESC
!
!     TSTRNG : input    time string
!
      CHARACTER (LEN=24) :: TSTRNG                                        40.02
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
      INTEGER    ITIME(6)
!
!  8. SUBROUTINE USED
!
!     DTTIST   (installation dependent subroutines)
!     DTINTI   (misc. routines)
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!                                                                         30.00
      CALL DTINTI (TIMESC, ITIME)
      CALL DTTIST (IOPT, TSTRNG, ITIME)
      DTTIWR = TSTRNG(1:18)                                               40.02
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE REPARM (NDSL, NDSD, IDLA, IDFM, RFORM,                   40.00
     &                   NHEDF, IDYN, NHEDT, LOGC, LOCAL, NHEDC)          40.95 40.00
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
!PUN      USE SIZES, ONLY: LOCALDIR                                           40.95
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
!     30.74: IJsbrand Haagsma (Include version)
!     30.80: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.00: Nico Booij (modifications)
!     40.02: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.10, Dec. 95: [fac] is read before 'fname' in view of later
!                     reading of tables; argument VFAC removed
!                     arguments LOGT, NHEDT, LOGC, NHEDC added
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.00, Jan. 98: SWAN specific statements modified; argument list
!                     changed
!     30.82, Sep. 98: Added type specification for HEDLIN
!     30.80, Dec. 98: Initialisation of NHEDT and NHEDC
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Sep. 00: Replaced computed GOTO by CASE construct
!     40.03, Jul. 00: TRIM used to improve readability of message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     reads parameters for reading an array from users input
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IDFM   : output   format index
!     IDLA   : output   lay-out indicator
!     IDYN   : input    indicate whether grid is dynamic or not
!     NDSD   : ??       unit number of the file from which to read the dataset
!     NDSL   : ??       unit number of the file containing the list of filenames
!     NHEDF  : output   number of heading lines in the file (once in each file)
!     NHEDT  : output   number of heading lines in the file before reading
!                       each time level
!     NHEDC  : output   number of heading lines in the file before each array
!                       or vector component
!
      INTEGER   IDFM, IDLA,  NDSL, NDSD, NHEDF, NHEDT, NHEDC, IDYN
!
!     LOCAL  : input    if True input file found in local directory
!     LOGC   : input    if True more than one component is read from file
!
      LOGICAL   LOGC, LOCAL
!
!     RFORM  : output   reading format
!
      CHARACTER RFORM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!
!     IENT   : Number of entries into this subroutine
!     IH     : ??
!     IOSTAT : input   0 : Full messages printed
!                      -1: Only error messages printed
!                      -2: No messages printed
!              output  error indicator
!
      INTEGER   IENT, IH, IOSTAT
!
!     HEDLIN : Content of a header line
!     KEYWIS : ??
!
      LOGICAL   KEYWIS, BNEW
!
!     OLDFIL : ??
!
      CHARACTER HEDLIN*80, OLDFIL *36
!                                                                         30.82
!  8. SUBROUTINE USED
!
      LOGICAL   STPNOW                                                    34.01
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE  IENT, OLDFIL
      DATA  IENT /0/
      DATA  OLDFIL /'                                    '/
      CALL STRACE (IENT, 'REPARM')
!
      CALL INKEYW ('STA', '   ')
      IF (KEYWIS('SERI')) THEN
        CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')
!       open namelist file and read first datafile name
        CALL FOR (NDSL, FILENM, 'OF', IOSTAT)                             40.00
        IF (STPNOW()) RETURN                                              34.01
        READ(NDSL, '(A36)') FILENM
      ELSE
        CALL INCSTR ('FNAME', FILENM, 'REQ', ' ')                         40.00
      ENDIF
!PUN      IF (LOCAL) FILENM = TRIM(LOCALDIR)//DIRCH2//TRIM(FILENM)            40.95
!
      IF (FILENM.NE.OLDFIL) THEN
        BNEW = .TRUE.
        NDSD = 0
        IDLA = 1
        IDFM = 0
        RFORM = ' '
        NHEDF = 0
        OLDFIL = FILENM
      ELSE
        BNEW = .FALSE.
      ENDIF
!
      CALL INKEYW ('STA', ' ')
      IF (BNEW) THEN
!       read lay-out indicator
        CALL ININTG ('IDLA', IDLA, 'UNC', 1)
!       names changed and order changed, ver 30.20 (Swan)
        CALL ININTG ('NHEDF', NHEDF, 'UNC', 0)                            30.20
        NHEDT = 0                                                         30.80
        IF (IDYN.GT.0) THEN
          CALL ININTG ('NHEDT', NHEDT, 'UNC', 0)                          30.21
        ENDIF
        NHEDC = 0                                                         30.80
        IF (LOGC) THEN
          CALL ININTG ('NHEDVEC', NHEDC, 'UNC', 0)                        30.20
        ENDIF
        CALL INKEYW ('STA', 'FREE')                                       30.06
        IDFM = 2
        IF (KEYWIS('FRE')) THEN
          IDFM = 0
        ELSE IF (KEYWIS('UNF')) THEN
          IDFM = -1
        ELSE IF (KEYWIS('FOR')) THEN
!         formatted read
          CALL ININTG ('IDFM', IDFM, 'NSKP', 2)
          SELECT CASE(IDFM)                                               40.02
          CASE(1)                                                         40.02
            RFORM = '(10X,12F5.0)'                                        40.02
          CASE(2)                                                         40.02
            CALL INCSTR ('FORM', RFORM, 'REQ', ' ')                       40.02
          CASE(5)                                                         40.02
            RFORM = '(16F5.0)'                                            40.02
          CASE(6)                                                         40.02
            RFORM = '(12F6.0)'                                            40.02
          CASE(8)                                                         40.02
            RFORM = '(10F8.0)'                                            40.02
          CASE DEFAULT                                                    40.02
            CALL MSGERR (2, 'illegal format number')                      40.02
            WRITE (PRINTF, 50) IDFM                                       40.02
  50        FORMAT (' -> ', I6)                                           40.02
          END SELECT                                                      40.02
        ELSE
          CALL WRNKEY                                                     30.06
          IDFM = 0
        ENDIF
!       --------------------------------------------------------
!                          open the file
!       --------------------------------------------------------
        IF (IDFM.NE.-1) THEN
          IOSTAT = 0
          CALL FOR (NDSD, FILENM, 'OF', IOSTAT)                           40.00
          IF (STPNOW()) RETURN                                            34.01
          IF (NHEDF.GT.0) THEN                                            40.00
            WRITE (PRINTF, '(A,A,A)') ' **  Heading lines file ',
     &      TRIM(FILENM), ' **'                                           40.03
            DO IH=1, NHEDF
              READ (NDSD, '(A80)') HEDLIN
              WRITE (PRINTF, '(A4,A80)') ' -> ', HEDLIN                   40.00
            ENDDO
          ENDIF
        ELSE
          IOSTAT = 0
          CALL FOR (NDSD, FILENM, 'OU', IOSTAT)                           40.00
          IF (STPNOW()) RETURN                                            34.01
          DO IH=1, NHEDF
            READ (NDSD)                                                   40.00
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE INAR2D (ARR, MXA, MYA, NDSL, NDSD, IDFM, RFORM,          40.00
     &  IDLA, VFAC, NHED, NHEDF)                                          40.00
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.08: Erick Rogers
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     01.05, Feb. 90: Before reading values in the array are divided by VFAC,
!                     in order to retain correct values for points where no
!                     value was given
!     01.06, Apr. 91: i/o status is printed if read error occurs
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.72, Sept 97: Corrected reading of heading lines for SERIES of files
!                     in dynamic mode
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     40.00, July 98: SWAN specific statements modified
!                     unformatted read: heading lines also read unformatted
!                     distinction between NDSD (data file) and NDSL (file list)
!     30.82, Sep. 98: Added INQUIRE statement to produce correct file name in
!                     case of a read error
!     34.01, Feb. 99: Introducing STPNOW
!     40.02, Sep. 00: Replaced computed GOTO with CASE construct
!     40.02, Sep. 00: Replaced reserved words IOSTAT with IOERR and STATUS with IERR
!     40.03, Jul. 00: END= added to READ statement for correct reading of series
!                     of files
!     40.03, Jul. 00: TRIM used to improve readability of message
!     40.13, Apr. 01: END=930 added in READ statement; corresponding error message added
!     40.08, Mar. 03: Changed an INQUIRE statement so that it does not produce
!                     misleading results.
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Reads a 2d array from dataset
!     is used to read e.g. bathymetry, one component of wind velocity
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IDFM   : input    format index
!     IDLAM  : input    lay-out indicator
!     MXA    : input    number of points along x-side of grid
!     MYA    : input    number of points along y-side of grid
!     NDSD   : input    unit number of the file from which to read the dataset
!     NDSL   : input    unit number of the file containing the list of filenames
!     NHEDF  : input    number of heading lines in the file (first lines).
!     NHEDL  : input    number of heading lines in the file
!                       before each array
!
      INTEGER   IDFM, IDLA, MXA, MYA, NDSD, NDSL, NHED, NHEDF
!
!     ARR    : input    results appear in this array
!     RFORM  : input    format used in reading data (char. string)
!     VFAC   : input    factor by which data must be multiplied.
!
      REAL      ARR(MXA,MYA), VFAC
!
      CHARACTER RFORM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IERR   : ??
!     IENT   : number of entries into this subroutine
!     IOERR  : input   0 : Full messages printed
!                      -1: Only error messages printed
!                      -2: No messages printed
!              output  error indicator
!     IH     : ??
!     IX     : ??
!     IY     : ??
!     NUMFIL : ??
!
      INTEGER   IERR, IENT, IOERR, IH, IX, IY, NUMFIL                     40.02
!
!     HEDLIN : Content of a header line
!
      CHARACTER HEDLIN *80
!
!  8. SUBROUTINE USED
!
      LOGICAL STPNOW                                                      34.01
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'INAR2D')
!
 999  IF (NDSD.LT.0) RETURN                                               40.00
!     no reading from file due to open error
!
!     *** NUMFIL is the number of that is open in one time step  **
      NUMFIL = 0                                                          30.00
      IF (ITEST.GE.100) THEN
        WRITE (PRINTF, 12) MXA, MYA, NDSD, IDFM, RFORM,                   40.00
     &  IDLA, VFAC, NHED
  12    FORMAT (' * TEST INAR2D *', 4I4, 1X, A16, I3, 1X, E12.4, I3)
      ENDIF
!
!     Read heading lines, and print the same:
!
  11  IF (NHED.GT.0) THEN
        IF (IDFM.LT.0) THEN                                               40.00
          IF (ITEST.GE.30)
     &             WRITE (PRINTF, '(I3,A)') NHED, ' Heading lines'        40.00
          DO 28 IH=1, NHED
            READ (NDSD, END=910)                                          40.03
  28      CONTINUE
        ELSE
          DO 30 IH=1, NHED
            READ (NDSD, '(A80)', end=910) HEDLIN                          40.03
            IF (IH.EQ.1) WRITE (PRINTF, '(A)') ' **  Heading lines  **'
            WRITE (PRINTF, '(A4,A80)') ' -> ', HEDLIN
  30      CONTINUE
        ENDIF
      ENDIF
!
!     divide existing values in the array by VFAC
!
      DO 39 IY = 1, MYA                                                   30.72
        DO 38 IX = 1, MXA
          ARR(IX,IY) = ARR(IX,IY) / VFAC
  38    CONTINUE                                                          30.72
  39  CONTINUE                                                            30.72
!
!     start reading of 2D-array
!
      IF (IDFM.EQ.0) THEN
!       free format read
        SELECT CASE(IDLA)                                                 40.02
        CASE(1)                                                           40.02
          DO IY=MYA, 1, -1                                                40.02
            READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                 40.02
     &      (ARR(IX,IY), IX=1,MXA)                                        40.02
          ENDDO                                                           40.02
        CASE(2)                                                           40.02
          READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                   40.02
     &               ((ARR(IX,IY), IX=1,MXA), IY=MYA,1,-1)                40.02
        CASE(3)                                                           40.02
          DO IY=1, MYA                                                    40.02
            READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                 40.02
     &               (ARR(IX,IY), IX=1,MXA)                               40.02
          ENDDO                                                           40.02
        CASE(4)                                                           40.02
          READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                   40.02
     &              ((ARR(IX,IY), IX=1,MXA), IY=1,MYA)                    40.02
        CASE(5)                                                           40.02
          DO IX=1, MXA                                                    40.02
            READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                 40.02
     &               (ARR(IX,IY), IY=1,MYA)                               40.02
          ENDDO                                                           40.02
        CASE(6)                                                           40.02
          READ (NDSD, *, END=910, ERR=920, IOSTAT=IERR)                   40.02
     &               ((ARR(IX,IY), IY=1,MYA), IX=1,MXA)                   40.02
        END SELECT                                                        40.02
      ELSE IF (IDFM.GT.0) THEN
!       read with fixed format
        SELECT CASE (IDLA)                                                40.02
        CASE(1)                                                           40.02
          DO IY=MYA, 1, -1                                                40.02
            READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)             40.02
     &      (ARR(IX,IY), IX=1,MXA)                                        40.02
          ENDDO                                                           40.02
        CASE(2)                                                           40.02
          READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)               40.02
     &    ((ARR(IX,IY), IX=1,MXA), IY=MYA,1,-1)                           40.02
        CASE(3)                                                           40.02
          DO IY=1, MYA                                                    40.02
            READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)             40.02
     &      (ARR(IX,IY), IX=1,MXA)                                        40.02
          ENDDO                                                           40.02
        CASE(4)                                                           40.02
          READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)               40.02
     &    ((ARR(IX,IY), IX=1,MXA), IY=1,MYA)                              40.02
        CASE(5)                                                           40.02
          DO IX=1, MXA                                                    40.02
            READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)             40.02
     &      (ARR(IX,IY), IY=1,MYA)                                        40.02
          ENDDO                                                           40.02
        CASE(6)                                                           40.02
          READ (NDSD, RFORM, END=910, ERR=920, IOSTAT=IERR)               40.02
     &    ((ARR(IX,IY), IY=1,MYA), IX=1,MXA)                              40.02
        END SELECT                                                        40.02
      ELSE
!       unformatted read
        SELECT CASE(IDLA)
        CASE(1)
          DO IY=MYA, 1, -1                                                40.02
            READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                    40.02
     &      (ARR(IX,IY), IX=1,MXA)                                        40.02
          ENDDO                                                           40.02
        CASE(2)                                                           40.02
          READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                      40.02
     &    ((ARR(IX,IY), IX=1,MXA), IY=MYA,1,-1)                           40.02
        CASE(3)                                                           40.02
          DO IY=1, MYA                                                    40.02
            READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                    40.02
     &      (ARR(IX,IY), IX=1,MXA)                                        40.02
          ENDDO                                                           40.02
        CASE(4)                                                           40.02
          READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                      40.02
     &    ((ARR(IX,IY), IX=1,MXA), IY=1,MYA)                              40.02
        CASE(5)                                                           40.02
          DO IX=1, MXA                                                    40.02
            READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                    40.02
     &      (ARR(IX,IY), IY=1,MYA)                                        40.02
          ENDDO                                                           40.02
        CASE(6)                                                           40.02
          READ (NDSD, END=910, ERR=920, IOSTAT=IERR)                      40.02
     &    ((ARR(IX,IY), IY=1,MYA), IX=1,MXA)                              40.02
        END SELECT                                                        40.02
      ENDIF
      GOTO 900                                                            40.02
!
!     *** End of data file, in case SERIES next file is opened
!     *** unit = NDSD is closed before the next one is opened
!
 910  CONTINUE
      CLOSE(NDSD)
      NUMFIL = NUMFIL + 1
      IF (NUMFIL .GE. 2) GO TO 911
      IF (NDSL.GT.0) THEN
        READ (NDSL, '(A)', END=930) FILENM                                40.13
        IF (IDFM.NE.-1) THEN
          IOERR = 0
          CALL FOR (NDSD, FILENM, 'OF', IOERR)                            40.02
          IF (STPNOW()) RETURN                                            34.01
        ELSE
          IOERR = 0
          CALL FOR (NDSD, FILENM, 'OU', IOERR)                            40.02
          IF (STPNOW()) RETURN                                            34.01
        ENDIF
!
!       Read heading lines, and print these:
!                                                                         30.72
  2     IF (NHEDF.GT.0) THEN                                              30.72
          IF (IDFM.LT.0) THEN                                             40.00
            IF (ITEST.GE.30) WRITE (PRINTF, '(I3,A,A)') NHEDF,
     &            ' Heading lines at begin of file ', TRIM(FILENM)        40.03
            DO 828 IH=1, NHEDF                                            40.00
              READ (NDSD)                                                 40.00
 828        CONTINUE
          ELSE                                                            40.00
            WRITE (PRINTF, '(A,A,A)') ' **  Heading lines file ',
     &      TRIM(FILENM), ' **'                                           40.03
            DO 830 IH=1, NHEDF                                            30.72
              READ (NDSD, '(A80)') HEDLIN                                 30.72
              WRITE (PRINTF, '(A4,A80)') ' -> ', HEDLIN                   30.72
 830        CONTINUE                                                      30.72
          ENDIF                                                           40.00
        ENDIF                                                             30.72
        GO TO 11
      ENDIF
!
!     error message when end of file is encountered
!
!     --- initialize FILENM so that previous value is not used            40.08
!         in case unit NDSD does not exist                                40.08
 911  FILENM='DUMMY'
!     --------------------------------------------------------------------40.08
!     THIS INQUIRE STATEMENT IS PROBLEMATIC, SINCE (AT LEAST              40.08
!     SOMETIMES) NDSD HAS ALREADY BEEN CLOSED, SO THE INQUIRE             40.08
!     STATEMENT SHOULD NOT WORK.                                          40.08
!     --------------------------------------------------------------------40.08
      INQUIRE (UNIT=NDSD, NAME=FILENM)
      CALL MSGERR (2, 'Unexpected end of file while reading '//
     &                 TRIM(FILENM))                                      40.13
      NDSD = 0                                                            40.00
      IDLA = -1
!     Value of IDLA=-1 signals end of file to calling program
!
      GOTO 900
!
!     --- initialize FILENM                                               40.08
 920  FILENM='DUMMY'                                                      40.08
      INQUIRE (UNIT=NDSD, NAME=FILENM)                                    30.82 40.08
      CALL MSGERR (2, 'Error while reading file '//TRIM(FILENM))          40.13
      WRITE (PRINTF, 922) IERR                                            40.02
 922  FORMAT (' i/o status ', I6)                                         40.00
      IDLA = -2                                                           40.00
!     Value of IDLA=-2 signals read error to calling program
!
!     Multiply all values in the array by VFAC
!
 900  DO 909 IY = 1, MYA                                                  30.72
        DO 908 IX = 1, MXA
          ARR(IX,IY) = ARR(IX,IY) * VFAC
 908    CONTINUE                                                          30.72
 909  CONTINUE                                                            30.72
!
 990  IF (ITEST.GE.100 .OR. IDLA.LT.0) THEN
        DO 996 IY=MYA, 1, -1
          WRITE (PRINTF, 994) (ARR(IX,IY), IX=1,MXA)
 994      FORMAT ((1X, 10E12.4))
 996    CONTINUE
      ENDIF
      RETURN

!     No more files in NDSL:
!     --- initialize FILENM                                               40.08
 930  FILENM='DUMMY'                                                      40.08
      INQUIRE (UNIT=NDSL, NAME=FILENM)                                    40.13 40.08
      CALL MSGERR (2, 'Series of input files ended in '//TRIM(FILENM))    40.13
      RETURN                                                              40.13

      END subroutine INAR2D
!*****************************************************************
!                                                                *
      SUBROUTINE STRACE (IENT, SUBNAM)
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!  0. AUTHORS
!
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This subroutine produces depending on the value of 'ITRACE'
!     a message containing the name 'SUBNAM'. the purpose of this
!     action is to detect the entry of a subroutine.
!
!  3. METHOD
!
!     the first executable statement of subroutine 'AAA' has to
!     be : CALL STRACE(IENT,'AAA')
!     further is necessary : DATA IENT/0/
!     IF ITRACE=0, no message
!     IF ITRACE>0, a message is printed up to ITRACE times
!
!  4. ARGUMENT VARIABLES
!
!     IENT   :  i/o    Number of entries into the calling subroutine
!
      INTEGER IENT
!
!     SUBNAM :  inp    name of the calling subroutine.
!
      CHARACTER SUBNAM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!                                                                         40.31
!$    LOGICAL,EXTERNAL :: OMP_IN_PARALLEL                                 40.31
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      IF (ITRACE.EQ.0) RETURN
      IF (IENT.GT.ITRACE) RETURN
!$    IF (OMP_IN_PARALLEL()) THEN                                         40.31
!$OMP MASTER                                                              40.31
!$       IENT=IENT+1                                                      40.31
!$       WRITE (PRTEST, 10) SUBNAM                                        40.31
!$       IF (SCREEN.NE.PRINTF) WRITE (SCREEN, 10) SUBNAM                  40.31
!$OMP END MASTER                                                          40.31
!$    ELSE                                                                40.31
         IENT=IENT+1
         WRITE (PRTEST, 10) SUBNAM
         IF ( SCREEN.NE.PRINTF .AND. IAMMASTER )                          40.95 40.30
     &                                         WRITE (SCREEN, 10) SUBNAM  40.30
!$    ENDIF                                                               40.31
  10  FORMAT (' ++ trace subr: ',A)
      RETURN
!  *  END OF SUBR. STRACE  *
      END
!*****************************************************************
!                                                                *
      SUBROUTINE MSGERR (LEV,STRING)
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
!  0. AUTHORS
!
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.03, Aug. 00: variable ERRFNM introduced in order to get correct
!                     message on UNIX system
!     40.02, Sep. 00: Removed STOP statement
!     40.13, Nov. 01: OPEN statement instead of CALL FOR
!                     to prevent recursive subroutines calling
!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Error messages are produced by subroutine MSGERR. if necessary
!     the value of LEVERR is increased.
!     In case of a high error level an error message file is opened
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     LEV    : indicates how severe the present error is
!     STRING : contents of the present error message
!
      INTEGER   LEV
!
      CHARACTER STRING*(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IERR   : if non-zero error message file was already opened unsuccessfully
!     IERRF  : unit reference number of the error message file
!     ILPOS  : actual length of error message filename
!
      INTEGER, SAVE :: IERR=0, IERRF=0                                    40.03
      INTEGER ILPOS                                                       40.30
!
!     ERRM   : error message prefix
!
      CHARACTER (LEN=17) :: ERRM                                          40.03
!
!     ERRFNM : name of error message file
!
      CHARACTER (LEN=LENFNM), SAVE :: ERRFNM = 'Errfile'                  40.31 40.03
!
!  8. SUBROUTINE USED
!
!     ---
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
!
      IF (LEV.GT.LEVERR) LEVERR=LEV
      IF (LEV.EQ.0) THEN
        ERRM = 'Message          '
      ELSE IF (LEV.EQ.1) THEN
        ERRM = 'Warning          '
      ELSE IF (LEV.EQ.2) THEN
        ERRM = 'Error            '
      ELSE IF (LEV.EQ.3) THEN
        ERRM = 'Severe error     '
      ELSE
        ERRM = 'Terminating error'
      ENDIF
      WRITE (PRINTF,12) ERRM, STRING
  12  FORMAT (' ** ', A, ': ',A)
      IF (LEV.GT.MAXERR) THEN
        IF (IERRF.EQ.0) THEN
          IF (IERR.NE.0) RETURN
!
!         append node number to ERRFNM in case of                         40.30
!         parallel computing                                              40.30
!
!          IF (PARLL) THEN                                                 40.30
!             ILPOS = INDEX ( ERRFNM, ' ' )-1                              40.30
!             WRITE(ERRFNM(ILPOS+1:ILPOS+4),13) INODE                      40.30
!  13         FORMAT('-',I3.3)                                             40.30
!          END IF                                                          40.30
!
!PUN          ERRFNM = TRIM(LOCALDIR)//DIRCH2//TRIM(ERRFNM)                   40.95
!PUN!
!          IERRF = 17                                                      40.13
!          OPEN (UNIT=IERRF, FILE=ERRFNM, FORM='FORMATTED')                40.13
        ENDIF
!        WRITE (IERRF,14) ERRM, STRING
!  14    FORMAT (A, ': ',A)
      ENDIF
!
      RETURN
!
      END SUBROUTINE MSGERR
!
!*****************************************************************
!                                                                *
      LOGICAL FUNCTION STPNOW()                                           30.82
!                                                                *
!*****************************************************************
!
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
!     30.82, Feb. 99: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.82: New function
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Function determines wheter the SWAN program should be stopped
!     due to a terminating error
!
!  3. Method
!
!     Compares two common variables (the maximum allowable error-level,
!     MAXERR and the actual error-level: LEVERR).
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT  : Number of entries into this subroutine
!
      INTEGER IENT
!
!  8. SUBROUTINE USED
!
!$    LOGICAL,EXTERNAL :: OMP_IN_PARALLEL
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'STPNOW')
!
      IF (LEVERR .GE. 4) THEN
        STPNOW = .TRUE.
      ELSE
        STPNOW = .FALSE.
      END IF
      IF (MAXERR.EQ.-1) STPNOW = .FALSE.
!$    IF (OMP_IN_PARALLEL()) STPNOW = .FALSE.
!
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE TABHED (PROGNM, LPR)
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!  0. AUTHORS
!
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.13, Jan. 01: VERTXT replaces VERNUM
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     prints the table heading, containing:
!     run description, 3 lines
!     name of institute, program name,
!     project name, run id.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     LPR    : input    unit ref. nr. for output
!
      INTEGER LPR
!
!     PROGNM : input    program name
!
      CHARACTER PROGNM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!  8. SUBROUTINE USED
!
!  9. SUBROUTINES CALLING
!
! 10. ERROR MESSAGES
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      WRITE (LPR, 10) PROJT1, INST
      WRITE (LPR, 20) PROJT2, PROGNM, VERTXT                              40.13
      WRITE (LPR, 30) PROJT3, PROJID, PROJNR
      WRITE (LPR, 40)
  10  FORMAT ('1', A72, ' | ', A40)
  20  FORMAT (1X,  A72, ' | ', A, '  version: ', A)                       40.13
  30  FORMAT (1X,  A72, ' | ', A16, 1X, A4)
  40  FORMAT (' --------------------------------------------------',
     &  '---------------------------------------------------------')
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE FOR (IUNIT, DDNAME, SF, IOSTAT)
!                                                                *
!*****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!     30.13: Nico Booij
!     30.70: Nico Booij
!     30.82: IJsbrand Haagsma
!     34.01: IJsbrand Haagsma
!     40.00, 40.03: Nico Booij
!     40.41: Marcel Zijlema
!     41.20: Casey Dietrich
!
!  1. Updates
!
!     30.13, Jan. 96: new structure
!     30.70, Feb. 98: terminating error if input file does not exist
!     30.82, Nov. 98: Introduced recordlength of 1000 for new files to
!                     avoid errors on the Cray-J90
!     34.01, Feb. 99: STOP statement removed
!     40.00, Feb. 99: DIRCH2 replaces DIRCH1 in filenames
!     40.03, May  00: modification for Linux: local copy of filename
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.20, Mar. 10: extension to tightly coupled ADCIRC+SWAN model
!
!  1. PURPOSE
!
!     General open file routine.
!
!  2. METHOD
!
!     FORTRAN 77 OPEN option.
!                INQUIRE
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!       IUNIT   int     input   =0 : get free unit number
!                               >0 : fixed unit number
!                       output  allocated unit number
!       DDNAME  char    input   ddname/filename string (empty if IUNIT>0)
!       SF      char*2  input   file qualifiers
!                               1st char: O(ld),N(ew),S(cratch),U(nknown)
!                               2nd char: F(ormatted),U(nformatted)
!       IOSTAT  int     input   0 : Full messages printed
!                               -1: Only error messages printed
!                               -2: No messages printed
!                       output  error indicator
!
      INTEGER   IUNIT, IOSTAT
      CHARACTER DDNAME*(LENFNM), SF*2                                     40.03
!
!  5. PARAMETER VAR. (CONSTANTS)
!
!     Error codes:
!
!       IOSTAT = IESUCC No errors
!       IOSTAT > 0      I/O error
!       IOSTAT = IENUNF No free unit number found
!       IOSTAT = IEUNBD Specified unit number out of bounds
!       IOSTAT = IENODD No filename supplied with IUNIT=0
!       IOSTAT = IEDDNM Incorrect filename supplied with IUNIT>0
!       IOSTAT = IEEXST Specified unit number does not exist
!       IOSTAT = IEOPEN Specified unit number already opened
!       IOSTAT = IESTAT Error in file qualifiers
!       IOSTAT = IENSCR Named scratch file
!       IOSTAT = IENSIO No specified I/O error
!
      INTEGER  IESUCC, IENUNF, IEUNBD, IENODD,
     &         IEDDNM, IEEXST, IEOPEN, IESTAT, IENSCR
      PARAMETER (IESUCC=  0,IENUNF= -1,IEUNBD= -2,IENODD= -3,
     &           IEDDNM= -4,IEEXST= -5,IEOPEN= -6,IESTAT= -7,
     &           IENSCR=-12)
!
!     EMPTY    blank string
!
      CHARACTER  EMPTY*(*)
      PARAMETER (EMPTY= '        ')
!
!  6. LOCAL VARIABLES
!
!     IENT      number of entries into this subroutine
!     IFO       format index
!     IFUN      free unit number
!     II        counter
!     IOSTTM    aux. error index
!     IS        file status index
!     IUTTM     aux. unit number
!
      INTEGER   IENT, IFO, IFUN, II, IOSTTM, IS, IUTTM
!
!     EXIST     if true, file exists
!     OPENED    if true, file is opened
!
      LOGICAL   EXIST, OPENED
!
!     S
!     F
!     FILTTM   auxiliary
!     FISTAT   file status, values: OLD, NEW, UNKNOWN
!     FORM     formatting, values: FORMATTED, UNFORMATTED
!     DDNAME_L local copy of DDNAME                                       40.03
!
      CHARACTER S, F, FILTTM *(LENFNM), DDNAME_L *(LENFNM)                40.03
      CHARACTER *11 FISTAT(4),FORM(2)
!
!  4. SUBROUTINES USED
!
!
!  5. ERROR MESSAGES
!
!       and error messages added using MSGERR
!
!
!  6. REMARKS
!
!       Free unit number search interval: FUNLO<=IUNIT<=FUNHI
!       FUNLO, FUNHI, IUNMIN and IUNMAX were initialized by OCPINI,
!       they are transmitted via module OCPCOMM4
!
!  7. STRUCTURE
!
!       ----------------------------------------------------------------
!       Check file qualifiers
!       ----------------------------------------------------------------
!       If IUNIT = 0
!       Then If DDNAME = ' '
!            Then error message
!            Else Inquire to find if file exists and is opened,
!                 and if so, to find correct unit number
!                 If file is not opened
!                 Then get a free unit number, assign value to IUNIT
!                      open the file
!                 Else assign correct unit number to IUNIT
!       Else Inquire to find if file exists and is opened,
!                   and if so, to find correct filename
!            If file with unit nr IUNIT is already open
!            Then If filename does not correspond to DDNAME
!                 Then Close file with old filename and unit IUNIT
!                      Open file with new filename DDNAME and unit IUNIT
!            Else If DDNAME is not empty
!                 Then Open file with new filename DDNAME and unit IUNIT
!                 Else Open file with unit IUNIT
!       ----------------------------------------------------------------
!
!  8. SOURCE TEXT
!
      SAVE      IENT, IFUN
!
      DATA FISTAT(1),FISTAT(2) / 'OLD','NEW'/
     &     FISTAT(3),FISTAT(4) / 'SCRATCH','UNKNOWN'/
     &     FORM(1),FORM(2) / 'FORMATTED','UNFORMATTED'/
!
      DATA IENT /0/, IFUN /0/
      CALL STRACE (IENT, 'FOR')
!
      IF (ITEST.GE.80) WRITE (PRTEST, 2) IUNIT, DDNAME, SF, IOSTAT
   2  FORMAT (' Entry FOR: ', I3, 1X, A36, A2, I7)
      DDNAME_L = DDNAME                                                   40.03
!
!     check file qualifiers
!
      IF ((IUNIT.NE.0) .AND.
     &   ((IUNIT .LT. IUNMIN) .OR. (IUNIT .GT. IUNMAX))) THEN
        IF (IOSTAT.GT.-2) CALL MSGERR (3, 'Unit number out of range')
        IOSTAT= IEUNBD
        RETURN
      END IF
!
      S   = SF(1:1)
      F   = SF(2:2)
      IS  = INDEX('ONSU',S)
      IFO = INDEX('FU',F)
      IF ((IS .EQ. 0) .OR. (IFO .EQ. 0)) THEN
        IF (IOSTAT.GT.-2) CALL MSGERR (3,'Error in file qualifiers')
        IOSTAT= IESTAT
        RETURN
      END IF
!
      IF ((S.EQ.'S').AND.(DDNAME.NE.EMPTY)) THEN
        IF (IOSTAT.GT.-2) CALL MSGERR (3, 'Named scratch file')
        IOSTAT= IENSCR
        RETURN
      END IF
!
      IF (DDNAME.NE.EMPTY) THEN                                           40.00
!       directory separation character is replaced in filenames           40.00
        DO II = 1, LEN(DDNAME)
          IF (DDNAME(II:II).EQ.DIRCH1) DDNAME(II:II) = DIRCH2             40.00
        ENDDO
      ENDIF
!
      IF (IUNIT .EQ. 0) THEN
         IF (DDNAME.EQ.EMPTY) THEN
            IF (IOSTAT.GT.-1) CALL MSGERR (3, 'No filename given')
            IOSTAT= IENODD
            RETURN
         ELSE
!           Was the file opened already ?
            INQUIRE (FILE=DDNAME, IOSTAT=IOSTTM, EXIST=EXIST,
     &      OPENED=OPENED, NUMBER=IUTTM)
            IF (IOSTTM .NE. IESUCC) THEN
               IF (IOSTAT.GT.-1) CALL MSGERR (2,
     &               'Inquire failed, filename: '//DDNAME_L)              40.03
               IOSTAT = IOSTTM
               RETURN
            ENDIF
!           If file does not exist, print term. error
            IF (IS.EQ.1 .AND. .NOT. EXIST) THEN                           30.70
               CALL MSGERR (4,
     &         'File cannot be opened/does not exist: '//DDNAME_L)        40.03
               IOSTAT = IEEXST
               RETURN
            END IF
            IF (OPENED) THEN
               IF (IOSTAT.GT.-1)
     &         CALL MSGERR (2, 'File is already opened: '//DDNAME_L)      40.03
               IOSTAT = IEOPEN
               IUNIT = IUTTM
               RETURN
            ENDIF
!ADC  60        CONTINUE                                                      41.20
!           Assign free unit number
            IF (IFUN.EQ.0) THEN
               IFUN = FUNLO
            ELSE
               IFUN = IFUN + 1
            ENDIF
!ADC!           the following unit numbers are reserved for ADCIRC            41.20
!ADC            IF ( IFUN.EQ.13 ) GOTO 60                                     41.20
!ADC            IF ( IFUN.EQ.14 ) GOTO 60                                     41.20
!ADC            IF ( IFUN.EQ.15 ) GOTO 60                                     41.20
!ADC            IF ( IFUN.EQ.22 ) GOTO 60                                     41.20
            IUNIT = IFUN
            IF (IUNIT .GT. FUNHI) THEN
               IF (IOSTAT.GT.-2) CALL MSGERR (3, 'All free units used')
               IOSTAT= IENUNF
            ENDIF
         END IF
         OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,FILE=DDNAME,              30.82
!/Cray     &         RECL=1000,                                                 30.82
!/SGI     &         RECL=1000,                                                 30.82
!CVIS     &         SHARED,                                                    40.41
     &         STATUS=FISTAT(IS),ACCESS='SEQUENTIAL',FORM=FORM(IFO))      30.82
      ELSE
         INQUIRE (UNIT=IUNIT, NAME=FILTTM, IOSTAT=IOSTTM,
     &            EXIST=EXIST, OPENED=OPENED)
         IF (IOSTTM .NE. IESUCC) THEN
            IF (IOSTAT.GT.-1) CALL MSGERR (2,
     &            'Inquire failed, filename: '//FILTTM)
            IOSTAT = IOSTTM
            RETURN
         ENDIF
         IF (OPENED) THEN
            IF (IOSTAT.GT.-1) THEN
              CALL MSGERR (1,
     &                   'File is already opened, filename: '//FILTTM)
            ENDIF
            IF (FILTTM.NE.DDNAME .AND. FILTTM.NE.EMPTY) THEN
              IF (IOSTAT.GT.-2) THEN
                WRITE (PRINTF, '(A, I4, 6A)') ' unit', IUNIT,
     &                 ' filenames: ', FILTTM, ' and: ', DDNAME
                CALL MSGERR (2, 'filename and unit number inconsistent')
              ENDIF
              IOSTAT = IEDDNM
!             close old file and open new one with given filename
              CLOSE (IUNIT)
              OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),
!/Cray     &         RECL=1000,                                                 30.82
!/SGI     &         RECL=1000,                                                 30.82
!CVIS     &         SHARED,                                                    40.41
     &              FILE=DDNAME,ACCESS='SEQUENTIAL',FORM=FORM(IFO))
              IF (IOSTTM.NE.IESUCC) IOSTAT = IOSTTM
              GOTO 80
            ENDIF
            IOSTAT = IEOPEN
            RETURN
         END IF
         IF (DDNAME.NE.EMPTY) THEN
            OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),
!/Cray     &         RECL=1000,                                                 30.82
!/SGI     &         RECL=1000,                                                 30.82
!CVIS     &         SHARED,                                                    40.41
     &      FILE=DDNAME,ACCESS='SEQUENTIAL',FORM=FORM(IFO))
         ELSE
            OPEN (UNIT=IUNIT,ERR=999,IOSTAT=IOSTTM,STATUS=FISTAT(IS),
!/Cray     &         RECL=1000,                                                 30.82
!/SGI     &         RECL=1000,                                                 30.82
!CVIS     &         SHARED,                                                    40.41
     &      ACCESS='SEQUENTIAL',FORM=FORM(IFO))
         END IF
      END IF
      HIOPEN = IFUN
  80  IF (ITEST.GE.30) WRITE (PRINTF, 82) IUNIT, DDNAME, SF
  82  FORMAT (' File opened: ', I6, 2X, A36, 2X, A2)
      RETURN
!
!     in case file cannot be opened:
!
 999  IF (IOSTAT.GT.-2) THEN
        CALL MSGERR (3, 'File open failed, filename: '//DDNAME_L)         40.03
        WRITE (PRINTF,15) DDNAME, IOSTTM, SF
  15    FORMAT (' File -> ', A36, 2X, ' IOSTAT=', I6, 4X, A2)
      ENDIF
      IUNIT = -1
      IOSTAT= IOSTTM
      RETURN
!  *  end of subroutine  FOR  *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION EQREAL (REAL1, REAL2 )                             30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!     30.72 IJsbrand Haagsma
!     30.60 Nico Booij
!     40.04 Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Oct. 97: Changed from EXCYES to make floating point point comparisons
!     30.60, July 97: new subroutine (EXCYES)
!     40.04, Aug. 00: introduced EPSILON and TINY
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     to determine whether a value (usually a value read from file)
!     is an exception value or not
!     Later (30.72) used to make comparisons of floating points within reasonable bounds
!
!  3. Method (updated...)
!
!     Checks whether ABS(REAL1-REAL2) .LE. TINY(REAL1) or whether this    40.04
!     difference is .LE. then EPS (= EPSILON(REAL1)*ABS(REAL1-REAL2) )    40.04
!
!  4. Argument variables
!
!     REAL1  : input    value that is to be tested
!     REAL2  : input    given exception value
!
      REAL      REAL1, REAL2
!
!  5. Parameter variables
!
!  6. Local variables
!
!     EPS    : Small number (related to REAL1 and its difference with REAL2)
!     IENT   : Number of entries into this subroutine
!
      REAL      EPS
      INTEGER   IENT
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWREAD
!     SWDIM
!     SIRAY
!     SWBOUN
!     SWODDC
!     SWOEXD
!     SWOEXA
!     SWOEXF
!     SWPLOT
!     SWSPEC
!     ISOLIN
!     SNYPT2
!     INCTIM
!     INDBLE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'EQREAL')
      EQREAL = .FALSE.
!
      EPS = EPSILON(REAL1)*ABS(REAL1-REAL2)                               40.04
      IF (EPS ==0) EPS = TINY(REAL1)                                      40.04
      IF (ABS(REAL1-REAL2) .GT. TINY(REAL1)) THEN                         40.04
        IF (ABS(REAL1-REAL2) .LT. EPS) EQREAL = .TRUE.                    40.04
      ELSE                                                                40.04
        EQREAL = .TRUE.                                                   40.04
      ENDIF                                                               40.04
      RETURN
!     end of subroutine EQREAL
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE LSPLIT(RELINE, DATITM, NUMITM)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
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
!  0. AUTHORS
!
!     40.00, 40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.00, Jan. 98: New subroutine for SWAN
!     40.03, Jun. 00: declaration updated, TRIM added for readability
!                     test output added
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     a line read from a file is separated into single data items
!     each data item is found in a string DATITM
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     NUMITM : input    max number of data items in array
!
      INTEGER, INTENT(IN) :: NUMITM                                       40.03
!
!     DATITM : output   array of data items
!     RELINE : input    string (read from an input file)
!
      CHARACTER (LEN=*), INTENT(OUT) :: DATITM(NUMITM)                    40.03
      CHARACTER (LEN=*), INTENT(IN) ::  RELINE                            40.03
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     CRL    : a character of the input line RELINE
!     QUOTE  : ' i.e. string delimiter
!
      CHARACTER QUOTE *1, CRL *1
!
!     ICR1   : ??
!     IENT   : Number of entries into this subroutine
!     ILL    : sequence number of character being processed
!     IITM   : counter of data items
!     LENLIN : lenght of an input line
!     RITM   : type of data, 0: empty string, 2: string enclosed
!              in quotes, 1: other
!
      INTEGER   ICR1, IENT, ILL, IITM, LENLIN, RITM
!
!     LCHSTR : if True, program is reading a string (enclosed in quotes)
!
      LOGICAL   LCHSTR
!
!  8. SUBROUTINE USED
!
!     ------
!
!  9. SUBROUTINES CALLING
!
!     SWBOUN
!
! 10. ERROR MESSAGES
!
!     Too many data items on input line
!
! 11. REMARKS
!
! 12. STRUCTURE
!
! 13. SOURCE TEXT
!
      SAVE      IENT, QUOTE
      DATA      IENT /0/, QUOTE /''''/
      CALL STRACE (IENT, 'LSPLIT')
!
      LENLIN = LEN(RELINE)
      LCHSTR = .FALSE.
      RITM   = 1
      ICR1   = 1
      DO IITM = 1, NUMITM
        DATITM(IITM) = '    '
      ENDDO
      IF (ITEST.GE.150) WRITE (PRTEST,*) ' test LSPLIT ', RELINE
!
!     free format: separate the line into data items
!     blanks and commas serve as separation between data items
!     DATITM is string containing one data item
!
      IITM = 1
      DO 170 ILL = 1, LENLIN
         CRL = RELINE(ILL:ILL)
         IF (LCHSTR) THEN
!           reading a character string enclosed in quotes
            IF (CRL.EQ.QUOTE) THEN
!              closing quote
               LCHSTR = .FALSE.
               RITM   = 2
               IF (IITM.GT.NUMITM) THEN
                 CALL MSGERR (2, 'too many items on input line')
                 WRITE (PRINTF, *) ' -> ', TRIM(RELINE)
               ENDIF
               DATITM(IITM) = RELINE (ICR1:ILL-1)
            ENDIF
         ELSE
            IF (CRL.EQ.',') THEN
               IF (RITM.EQ.0) THEN
!                 empty item
                  IITM = IITM + 1
                  IF (IITM.GT.NUMITM) THEN
                    CALL MSGERR (2, 'too many items on input line')
                    WRITE (PRINTF, *) ' -> ', TRIM(RELINE)                40.03
                  ENDIF
                  DATITM(IITM) = '    '
               ELSE
                  IF (RITM.EQ.1) DATITM(IITM) = RELINE(ICR1:ILL)
                  RITM = 0
               ENDIF
            ELSE IF (CRL.EQ.' ' .OR. CRL.EQ.TABC) THEN
               IF (RITM.EQ.1) THEN
                  IF (IITM.GT.NUMITM) THEN
                    CALL MSGERR (2, 'too many items on input line')
                    WRITE (PRINTF, *) ' -> ', TRIM(RELINE)                40.03
                  ENDIF
                  DATITM(IITM) = RELINE(ICR1:ILL)
                  RITM = 2
               ENDIF
            ELSE
               IF (RITM.NE.1) THEN
                  IITM = IITM + 1
                  IF (IITM.GT.NUMITM) THEN
                    CALL MSGERR (2, 'too many items on input line')
                    WRITE (PRINTF, *) ' -> ', TRIM(RELINE)                40.03
                  ENDIF
                  IF (CRL.EQ.QUOTE) THEN
                     ICR1 = ILL+1
                     LCHSTR = .TRUE.
                  ELSE
                     ICR1 = ILL
                     RITM = 1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (ITEST.GE.250) WRITE (PRTEST, 165) CRL, RITM,
     &                IITM, ICR1
 165     FORMAT (' test LSPLIT ', A1, 3I3, 2X, A20)
 170  CONTINUE
      IF (ITEST.GE.130) THEN                                              40.03
        DO IITM = 1, NUMITM
          WRITE (PRTEST, 810) IITM, DATITM(IITM)
 810      FORMAT (' LSPLIT data item ', I2, ' is: ', A)
        ENDDO
      ENDIF
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE BUGFIX (FIXABC)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM2                                                        40.41
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
!     40.03  Nico Booij
!
!  1. UPDATE
!
!     40.03, May  00: new subroutine
!
!  2. PURPOSE
!
!     Adding one character to the version character string
!
!  3. METHOD
!
!
!  4. Argument variables
!
!       FIXABC  char   input    character indicating a bugfix

      CHARACTER (LEN=1), INTENT(IN) :: FIXABC
!
!  5. Parameter variables
!
!  6. Local variables
!
!       IC      counter of characters
!
      INTEGER   IC
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       for characters in VERTXT starting at end, do
!           if character is not blank
!           then replace previous character by FIXABC
!       ----------------------------------------------------------------
!
! 13. Source text
!
      DO IC = LEN(VERTXT), 1, -1
        IF (VERTXT(IC:IC) .NE. ' ') THEN
          VERTXT(IC+1:IC+1) = FIXABC
          GOTO 80
        ENDIF
      ENDDO
  80  RETURN
!     end of subroutine BUGFIX
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE COPYCH (STRING, MOVE, IARRAY, LENARR, IERR)              30.81
!                                                                      *
!***********************************************************************
!
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
!  0. AUTHORS
!
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     ver 30.01
!     30.81, Nov. 98: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     30.81, Jan. 99: Replaced variable FROM by FROM_ and TO by TO_ (because
!                     FROM and TO are reserved words)
!     40.03, Nov. 99: LENS2 removed from WRITE statement (value not yet known)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     copy a string into an integer array or vice-versa
!     MOVE (TO_ or FROM_) indicates copying direction                     30.81
!
!  3. METHOD
!
!     ---
!
!  4. ARGUMENT VARIABLES
!
!     IARRAY : output   an integer array
!     LENARR : input    length of array IARRAY
!     IERR   : output   error status: 0=no error, 9=end-of-file           30.81
!
      INTEGER   IARRAY(*), LENARR, IENT,                                  30.72
     &          IERR                                                      30.81
!
!     STRING : i/o      a character string
!     MOVE   : input    if MOVE=TO_, STRING is copied to IARRAY           30.81
!                       if MOVE=FROM_, STRING is copied from IARRAY       30.81
!
      CHARACTER MOVE *1, STRING *(*)
!
!  5. PARAMETER VARIABLES
!
!     OPMLFC : largest allowed integer character (ASCII) code + 1
!     OPMNLI : number of characters that can be stored in one integer number
!
      INTEGER   OPMLFC, OPMNLI                                            30.72
!
      PARAMETER (OPMNLI=4, OPMLFC=128)
!
!  6. LOCAL VARIABLES
!
!     IC     : counter
!     IENT   : number of entries into this subroutine
!     II     : counter
!     LENS1  : length of a string
!     LENS2  : length of a string
!     LL     : integer representation of a character
!     MC1    : integer converted to/from character
!     MCHAR  : integer converted to/from character
!     MM     : aux. number
!     NSL    : position of character in string
!
      INTEGER   IC, II, LENS1, LENS2, LL, MC1, MCHAR, MM, NSL
!
!     CC     : a single character
!     CHAR   : intrinsic character function, translates integer to character
!     FROM_  : 'F'
!     TO_    : 'T'
!
      CHARACTER CC, TO_, FROM_                                            30.81
!
!  8. SUBROUTINE USED
!
!     CHAR, ICHAR (intrinsic functions)
!
!  9. SUBROUTINES CALLING
!
!     ---
!
! 10. ERROR MESSAGES
!
!     If PINDEX or PPLACE is out of range an error message is printed
!
! 11. REMARKS
!
!     ---
!
! 12. STRUCTURE
!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
! 13. SOURCE TEXT
!
      SAVE        IENT, TO_, FROM_                                        30.81
      DATA        IENT /0/, TO_ /'T'/, FROM_ /'F'/                        30.81
      CALL STRACE (IENT,'COPYCH')
!
      LENS1 = LEN(STRING)
      IF (LENARR.GT.360) THEN
        CALL MSGERR (2, 'extremely long string in COPYCH')
        WRITE (PRTEST, *) ' test COPYCH  ',
     &          MOVE, LENARR, LENS1, ' ', STRING(1:80)                    40.03
        LENARR = 360
      ENDIF
      LENS2 = LENARR*OPMNLI
!
      IF (MOVE .EQ. TO_) THEN                                             30.81
        NSL = 0
        DO 60 II = 1, LENARR
          MCHAR = 0
          DO 40 IC = 1, OPMNLI
            NSL = NSL + 1
            IF (NSL .LE. LENS1) THEN
              CC = STRING(NSL:NSL)
            ELSE
              CC = ' '
            ENDIF
            LL = ICHAR(CC)
            IF (LL.GE.OPMLFC) THEN
              IERR = 803                                                  30.81
              WRITE (PRTEST, 33) CC
  33          FORMAT (' character cannot be copied: ', A1)
              LL = ICHAR ('?')
            ENDIF
            MCHAR = OPMLFC*MCHAR + LL
!            IF (ITEST.GE.250) WRITE (PRTEST, *) NSL, CC, LL
  40      CONTINUE
          IARRAY(II) = MCHAR
  60    CONTINUE
        IF (LENS1.GT.LENS2) THEN
          DO 70 II = LENS1+1, LENS2
            IF (STRING(II:II) .NE. ' ') THEN
              IERR = 801                                                  30.81
              CALL MSGERR(1, 'string longer than capacity of array')
              IF (ITEST.GE.50) WRITE (PRTEST, *) ' test COPYCH  ',
     &               MOVE, LENARR, LENS1, LENS2, ' ', STRING(1:80)
              GOTO 165
            ENDIF
  70      CONTINUE
        ENDIF
      ELSE IF (MOVE .EQ. FROM_) THEN                                      30.81
!
!       character string copied from an array
!
!       first the string is filled with blanks
        STRING = '    '
        NSL = 0
        DO 160 II = 1, LENARR
          MC1 = IARRAY(II)
          DO 140 IC = 1, OPMNLI
            MM  = OPMLFC ** (OPMNLI-IC)
            LL  = MC1 / MM
            NSL = NSL + 1
            IF (NSL .LE. LENS1) THEN
              STRING(NSL:NSL) = CHAR(LL)
            ELSE
              IF (CHAR(LL) .NE. ' ') THEN
                IF (IERR.NE.802)                                          30.81
     &          CALL MSGERR(1, 'string shorter than capacity of array')
                IF (ITEST.GE.50) WRITE (PRTEST, *) ' test COPYCH  ',
     &               MOVE, LENARR, LENS1, LENS2, ' ', STRING
                IERR = 802                                                30.81
                GOTO 165
              ENDIF
            ENDIF
            MC1 = MC1 - LL * MM
!           IF (ITEST.GE.250) WRITE (PRTEST, *) NSL, LL, STRING(NSL:NSL)
 140      CONTINUE
          IF (MC1.NE.0) WRITE (PRINTF, *) ' Error COPYCH'
 160    CONTINUE
      ELSE
        CALL MSGERR (2, 'error COPYCH, argument MOVE')
      ENDIF
 165  IF (ITEST.GE.230) WRITE (PRTEST, 167) LENS1, STRING, MOVE,
     &         (IARRAY(II), II=1,LENARR)
 167  FORMAT (' exit COPYCH ', I3, 1X, A20, 1X, A1, 4(1X,I12))
      RETURN
!*    end of subroutine COPYCH   **
      END
