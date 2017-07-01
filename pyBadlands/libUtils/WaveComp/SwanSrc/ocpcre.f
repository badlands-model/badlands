!               OCEAN PACK  command reading routines
!
!  Contents of this file:
!     RDINIT
!     NWLINE
!     INKEYW
!     INREAL
!     INDBLE
!     ININTG
!     INCSTR
!     INCTIM
!     ININTV
!     LEESEL
!     GETKAR
!     PUTKAR
!     UPCASE
!     KEYWIS
!     WRNKEY
!     IGNORE
!     RDHMS
!
!****************************************************************
!                                                               *
      SUBROUTINE RDINIT
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Initialises the command reading system
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'RDINIT')
      KAR = ';'
      KARNR = LINELN + 1                                                  40.00
      ELTYPE = 'USED'
      BLANK = '    '
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE NWLINE
!                                                               *
!****************************************************************
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
!     34.01: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     34.01, Feb. 99: Changed STOP statement in a MSGERR(4,'message')
!     40.03, Apr. 99: length of command lines changed from 80 to LINELN (=120)
!                     name of input file included in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Jumps to reading of the next input line,
!     if the end of the previous one is reached.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA  IENT/0/
      CALL STRACE (IENT,'NWLINE')
   5  IF ((ELTYPE.EQ.'USED').OR.(ELTYPE.EQ.'EOR')) CALL LEESEL
      IF (ELTYPE.EQ.'EOF') GOTO 90
      IF (ELTYPE.EQ.'KEY' .AND. KEYWRD.NE.'        ') GOTO 50
      IF (ELTYPE.EQ.'INT') GOTO 50
      IF (ELTYPE.EQ.'REAL') GOTO 50
      IF (ELTYPE.EQ.'CHAR') GOTO 50
      IF (KARNR.LE.LINELN) GOTO 50                                        40.03
!     The end of the previous line is reached, there are no more
!     unprocessed data items on that line.
!     Jump to new line can take place.
      WRITE (PRINTF,9) '    '
   9  FORMAT (A4)
      KARNR=0
      KAR=' '
      ELTYPE='USED'
      GOTO 5
  90  IF (ITEST.GE.10) THEN
        INQUIRE (UNIT=INPUTF, NAME=FILENM)                                40.03
        WRITE (PRINTF, *) ' end of input file '//FILENM                   40.03
      ENDIF
  50  RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE INKEYW (KONT, CSTA)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     ver 30.70, Jan. 1998: data type 'OTHR' is condidered
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     this subroutine reads a keyword.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     KONT   : action to be taken if no keyword is found in input:
!              'REQ' (required) error message
!              'STA' (standard) the value of csta is assigned to keywrd.
!
!     CSTA   : see above.
!
      CHARACTER CSTA *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENS   : length of default string (CSTA)
!
      INTEGER   IENT, LENS
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
      SAVE IENT
      DATA  IENT/0/
      CALL  STRACE ( IENT, 'INKEYW')
!
!     if necessary, a new data item is read.
!
      IF (ELTYPE.EQ.'KEY' .AND. KEYWRD.EQ.'        ') GOTO 510
      IF (ELTYPE.EQ.'KEY') GOTO 900
      IF (ELTYPE.EQ.'EOR') GOTO 510
      IF (ELTYPE.EQ.'USED') GOTO 510
      GOTO 520
 510  CALL LEESEL
 520  IF (ELTYPE.EQ.'KEY') GOTO 900
!     KEYWORD IS READ
      IF ((KONT.EQ.'STA').OR.(KONT.EQ.'NSKP')) THEN
        LENS = LEN(CSTA)
        IF (LENS.GE.8) THEN
          KEYWRD = CSTA(1:8)
        ELSE
          KEYWRD = '        '
          KEYWRD(1:LENS) = CSTA
        ENDIF
        GOTO 900
      ENDIF
!     at the end of the input 'STOP' is generated.
      IF (ELTYPE.EQ.'EOF') THEN
        KEYWRD='STOP'
        CALL MSGERR (2, 'STOP statement is missing')
        GOTO 900
      ENDIF
!     ----------------------------------------------------------
!     Data appear where a keyword is expected.
!     The user must be informed.
!     ----------------------------------------------------------
      IF (ELTYPE.EQ.'EOR') THEN
        KEYWRD = '        '
        GOTO 900
      ENDIF
      IF (ELTYPE.EQ.'INT') THEN
        CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
        GOTO 510
      ENDIF
      IF (ELTYPE.EQ.'REAL') THEN
        CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
        GOTO 510
      ENDIF
      IF (ELTYPE.EQ.'CHAR' .OR. ELTYPE.EQ.'OTHR') THEN                    30.70
        CALL MSGERR (2, 'Data field skipped:'//ELTEXT)
        GOTO 510
      ENDIF
      IF (ELTYPE.EQ.'EMPT') THEN
        CALL MSGERR (2, 'Empty data field skipped')
        GOTO 510
      ENDIF
      CALL MSGERR (3, 'Error subr. INKEYW')
!     ----------------------------------------------------------
 900  IF (ITEST.GE.10) WRITE (PRINTF,910) KEYWRD
 910  FORMAT (' KEYWORD: ',A8)
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE INREAL (NAAM, R, KONT, RSTA)
!                                                               *
!****************************************************************
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
!     30.82: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     20.04, Aug. 93: logical CHGVAL is introduced it  is made True if
!                     user changes value of an input parameter via INREAL
!     30.82, Sep. 98: To avoid errors using the Cray-cf90 compiler
!                     introduced a dummy
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Reads a REAL number in free format.
!
!  3. METHOD
!
!     Uses the function INDBLE to read the number
!
!  4. ARGUMENT VARIABLES
!
!     R      : The value of the variable that is to be read.
!     RSTA   : Reference value needed for KONT='STA'or 'RQI'
!
      REAL      R, RSTA
!
!     KONT   : What to do with the variable?
!              ='REQ' : variable is required
!              ='UNC' : if no variable, then variable will not be changed
!              ='STA' : if no variable, then variable will get value of RSTA
!              ='RQI' : variable may not have the value of RSTA
!              ='REP' : (REPEAT)
!              ='NSKP': (NO SKIP) if data item is of different type,
!                       value is left unchanged
!     NAAM   : Name of the variable according to the user manual.
!
      CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!     DRSTA  : Double precision variant of RSTA
!     RDBL   : Double precision variant of R
!
      DOUBLE    PRECISION RDBL, DRSTA
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE ( IENT, 'INREAL')
!
      RDBL = DBLE(R)
      DRSTA  = DBLE(RSTA)                                                 30.82
      CALL INDBLE (NAAM, RDBL, KONT, DRSTA)                               30.82
!
!     RDBL may have changed due to the value of KONT
!
      R = REAL(RDBL)
      RETURN
!
!     End of subroutine INREAL
!
      END
!
!****************************************************************
!                                                               *
      SUBROUTINE INDBLE (NAAM, R, KONT, RSTA)
!                                                               *
!****************************************************************
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
!     30.72: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     30.72, Oct. 97: Introduced logical function EQREAL for floating point
!                     comparisons
!     20.05, Aug. 93: NEW subroutine for double prec. data
!     40.03, Feb. 00: local copy of NAAM used in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Reads a DOUBLE PRECISION number, in free format.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     R      : THE VARIABLE THAT IS TO BE READ.
!     RSTA   : SEE ABOVE
!
      DOUBLE    PRECISION R, RSTA
!
!     KONT   : What to do with the variable?
!              ='REQ'; Value in input file is required
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get value of RSTA
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!     NAAM   : name of the variable according to the user's manual.
!
      CHARACTER KONT *(*), NAAM *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENNM  : length of string NAAM
!
      INTEGER   IENT, LENNM
!
!     NAAM_L : local copy of NAAM                                         40.03
!
      CHARACTER (LEN=40) :: NAAM_L                                        40.03
!
!     EQREAL :
!
      LOGICAL   EQREAL                                                    30.72
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE ( IENT, 'INDBLE')
!
!     if necessary, a new data item is read
!
      CHGVAL = .FALSE.
      NAAM_L = NAAM                                                       40.03
      IF (ELTYPE.EQ.'USED') CALL  LEESEL
!     consider type of data item
      IF (ELTYPE.EQ.'KEY') THEN
!       find out whether in input is written: NAAM=...
        LENNM = LEN(NAAM)
        IF (NAAM.NE.ELTEXT(1:LENNM)) GOTO 20
        ELTYPE = 'USED'
        CALL LEESEL
        IF (ELTYPE.EQ.'REAL') GOTO 12
        IF (ELTYPE.EQ.'INT') GOTO 14
        GOTO 50
      ENDIF
!     if type is 'REAL' or 'INT', its value is assigned to R
  12  IF (ELTYPE.EQ.'REAL') THEN
        R      = ELREAL
        ELTYPE = 'USED'
        IF (.NOT.EQREAL(REAL(R),REAL(RSTA))) CHGVAL = .TRUE.              30.72
        GOTO 80
      ENDIF
  14  IF (ELTYPE.EQ.'INT') THEN
        R      = DBLE(ELINT)
        ELTYPE = 'USED'
        IF  (.NOT.EQREAL(REAL(R),REAL(RSTA))) CHGVAL = .TRUE.             30.72
        GOTO 80
      ENDIF
      IF (ELTYPE.EQ.'EOR') THEN
!       find out whether end of repeat is reached,
!       if so make R=RSTA.
        IF (KONT.NE.'REP') GOTO 20
        ELTYPE='USED'
        GOTO 70
      ENDIF
      IF (ELTYPE.EQ.'EOF') GOTO 20
      IF (ELTYPE.EQ.'EMPT') THEN
        ELTYPE='USED'
        GOTO 20
      ENDIF
      IF (ELTYPE.EQ.'ERR') GOTO 50
      IF (ELTYPE.EQ.'CHAR' .OR. ELTYPE.EQ.'OTHR') THEN                    30.04
        IF (KONT.EQ.'NSKP') GOTO 70
        CALL MSGERR (3,
     &           'Wrong type of data for variable '//NAAM_L)              40.03
        WRITE (PRINTF,18) NAAM, ELTEXT(1:LENCST)
  18    FORMAT (' -> ',A, '  item=', A)
        ELTYPE='USED'
        GOTO 70                                                           20.01
      ENDIF
      CALL MSGERR (3, 'Error subr. INREAL')
      WRITE (PRINTF, '(1X,A,A)') ELTYPE, KONT
!     -------------------------------------------------------------
!     data item of different type is read, action according to KONT
!     -------------------------------------------------------------
  20  IF (KONT.EQ.'REQ') GOTO 30
      IF (KONT.EQ.'REP') GOTO 70
      IF (KONT.EQ.'RQI') GOTO 28
      IF (KONT.EQ.'STA') GOTO 70
      IF (KONT.EQ.'NSKP') GOTO 70
      IF (KONT.EQ.'UNC') GOTO 80
!
      CALL MSGERR (3, 'Error subr. INREAL')
      WRITE (PRINTF, '(1X,A,A)') ELTYPE, KONT
      GOTO 70
!
  28  IF  (.NOT.EQREAL(REAL(R),REAL(RSTA))) GOTO 80                       30.72
  30  CALL MSGERR (3, 'No value for variable '//NAAM_L)                   40.03
      WRITE (PRINTF,18) NAAM, ELTEXT(1:LENCST)
      GOTO 70
!
  50  CALL MSGERR (3, 'Read error with variable '//NAAM_L)                40.03
      WRITE (PRINTF,18) NAAM, ELTEXT(1:LENCST)
      ELTYPE='USED'
!
  70  R      = RSTA
  80  IF (ITEST.GE.10) WRITE (PRINTF, 85) NAAM, R
  85  FORMAT (1X,A8,'=',D12.4)
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE ININTG (NAAM, IV, KONT, ISTA)
!                                                               *
!****************************************************************
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
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     20.04, August 1993:  logical CHGVAL is introduced
!                          it is made True if user changes value
!                          of an input parameter via ININTG
!     40.03, Feb. 00: local copy of NAAM used in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Reads an integer number, in free format
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IV     :  integer variable which is to be assigned a value
!     ISTA   :  default value
!
      INTEGER   IV, ISTA
!
!     NAAM   :  name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
      CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!     PARAMETERS: SEE SUBR. INREAL
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENNM  : length of the string NAAM
!
      INTEGER   IENT, LENNM
!
!     NAAM_L : local copy of NAAM                                         40.03
!
      CHARACTER (LEN=40) :: NAAM_L                                        40.03
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
      SAVE IENT
      DATA  IENT /0/
      CALL  STRACE ( IENT, 'ININTG')
!
      CHGVAL = .FALSE.
      NAAM_L = NAAM                                                       40.03
!     IF NECESSARY, A NEW DATA ITEM IS READ
      IF (ELTYPE.EQ.'USED') CALL  LEESEL
!     CONSIDER TYPE OF DAT ITEM.
      IF (ELTYPE.EQ.'EOR') GOTO 14
      IF (ELTYPE.EQ.'KEY') GOTO 16
      IF (ELTYPE.EQ.'EOF') GOTO 20
      IF (ELTYPE.EQ.'EMPT') GOTO 18
      IF (ELTYPE.EQ.'ERR') GOTO 50
      IF (ELTYPE.EQ.'CHAR') GOTO 40
      IF (ELTYPE.EQ.'OTHR') GOTO 40
      IF (ELTYPE.EQ.'REAL') GOTO 40
!     TYPE IS 'INT', VALUE IS ASSIGNED
  12  IV=ELINT
      IF (IV.NE.ISTA) CHGVAL = .TRUE.
      ELTYPE='USED'
      GOTO 80
  14  IF (KONT.NE.'REP') GOTO 20
      ELTYPE='USED'
      GOTO 70                                                             20.01
  16  LENNM = LEN(NAAM)
      IF (NAAM.NE.ELTEXT(1:LENNM)) GOTO 20
      CALL LEESEL
      IF (ELTYPE.EQ.'INT') GOTO 12
      GOTO 50
  18  ELTYPE='USED'
  20  IF (KONT.EQ.'REQ') GOTO 30
      IF (KONT.EQ.'REP') GOTO 70
      IF (KONT.EQ.'RQI') GOTO 28
      IF (KONT.EQ.'STA') GOTO 70
      IF (KONT.EQ.'NSKP') GOTO 70
      GOTO 80
  28  IF (IV.NE.ISTA) GOTO 80
  30  CALL MSGERR (2, 'No value for variable '//NAAM_L)                   40.03
      GOTO 70                                                             20.01
  40  IF (KONT.EQ.'NSKP') GOTO 70
      CALL MSGERR (2, 'Wrong type of data for variable '//NAAM_L)         40.03
      WRITE (PRINTF,41) NAAM, ELTEXT(1:LENCST)                            30.04
  41  FORMAT (' -> ',A8, '  item read=', A)
      ELTYPE='USED'
      GOTO 70                                                             20.01
  50  CALL MSGERR (2, 'Read error with variable '//NAAM_L)                40.03
      WRITE (PRINTF,41) NAAM, ELTEXT(1:LENCST)                            30.04
      ELTYPE='USED'
  70  IV=ISTA                                                             20.01
  80  IF (ITEST.GE.10) WRITE (PRINTF, 85) NAAM, IV
  85  FORMAT (1X,A8,'=',I6)
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE INCSTR (NAAM, C, KONT, CSTA)
!                                                               *
!****************************************************************
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
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     20.04, August 1993:  logical CHGVAL is introduced
!                          it is made True if user changes value
!                          of an input parameter via INCSTR
!     40.03, Feb. 00: local copy of NAAM used in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Reads a string in free format
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of CSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!     C      : string that is to be read from input file
!     CSTA   : default value of the string
!
      CHARACTER NAAM *(*), KONT *(*), C *(*), CSTA *(*)
!
!  5. PARAMETER VARIABLES
!
!     Parameters: see program documentation.
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENNM  : length of the string NAAM
!     LENW   : length of the string C
!     NS     : length of the string CSTA
!
      INTEGER   IENT, LENNM, LENW, NS
!
!     NAAM_L : local copy of NAAM                                         40.03
!
      CHARACTER (LEN=40) :: NAAM_L                                        40.03
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
      SAVE IENT
      DATA  IENT /0/
      CALL  STRACE ( IENT, 'INCSTR')
!
      CHGVAL = .FALSE.
      NAAM_L = NAAM                                                       40.03
      LENW = LEN(C)
!     IF NECESSARY, A NEW DATA ITEM IS READ.
      IF (ELTYPE.EQ.'USED') CALL  LEESEL
!     CONSIDER TYPE OF DATA ITEM.
      IF (ELTYPE.EQ.'KEY') THEN
!       FIND OUT WHETHER IN INPUT IS WRITTEN:  NAAM=....
        LENNM = LEN(NAAM)
        IF (NAAM.NE.ELTEXT(1:LENNM)) GOTO 20
        ELTYPE = 'USED'
        CALL LEESEL
        IF (ELTYPE.EQ.'CHAR') GOTO 12
        GOTO 50
      ENDIF
  12  IF (ELTYPE.EQ.'CHAR') THEN
!       TYPE IS 'CHAR', VALUE IS ASSIGNED.
        IF (LENCST.GT.LENW) THEN
          CALL MSGERR (2,
     &       'too long string given for: '//NAAM_L)                       40.03
          WRITE (PRINTF, 13) NAAM, ELTEXT(1:LENCST)
  13      FORMAT (' name=', A, ' string=', A)
        ENDIF
        C = ELTEXT(1:LENW)
        IF (C.NE.CSTA) CHGVAL = .TRUE.
        ELTYPE='USED'
        GOTO 80
      ENDIF
      IF (ELTYPE.EQ.'EOR') THEN
!       END OF REPEAT GROUP, STANDARD VALUE IS ASSIGNED.
        IF (KONT.NE.'REP') GOTO 20
        ELTYPE='USED'
        GOTO 70
      ENDIF
      IF (ELTYPE.EQ.'EOF') GOTO 20
      IF (ELTYPE.EQ.'EMPT') THEN
        ELTYPE='USED'
        GOTO 20
      ENDIF
      IF (ELTYPE.EQ.'ERR') GOTO 50
      IF (ELTYPE.EQ.'INT') GOTO 40
      IF (ELTYPE.EQ.'REAL') GOTO 40
      IF (ELTYPE.EQ.'OTHR') GOTO 40                                       30.04
      CALL MSGERR (3, 'Error subr. INCSTR')
      WRITE (PRINTF, '(1X,A,1X,A)') ELTYPE, KONT
      GOTO 80
!     --------------------------------------------------------
!     No value provided, action is taken according to 'KONT'.
!     --------------------------------------------------------
  20  IF (KONT.EQ.'REQ') GOTO 30
      IF (KONT.EQ.'REP') GOTO 70
      IF (KONT.EQ.'RQI') GOTO 28
      IF (KONT.EQ.'STA') GOTO 70
      IF (KONT.EQ.'NSKP') GOTO 70
      IF (KONT.EQ.'UNC') GOTO 80
      CALL MSGERR (3, 'Error subr. INCSTR')
      WRITE (PRINTF, '(1X,A,1X,A)') ELTYPE, KONT
      GOTO 70                                                             20.01
!
  28  IF (C(1:LENW).NE.CSTA(1:LENW)) GOTO 80
  30  CALL MSGERR (3, 'No value for variable '//NAAM_L)                   40.03
      GOTO 70                                                             20.01
  40  IF (KONT.EQ.'NSKP') GOTO 70
      CALL MSGERR (3, 'Wrong type of data for variable '//NAAM_L)         40.03
      WRITE (PRINTF,41) NAAM, ELTEXT(1:LENCST)
  41  FORMAT (' -> ',A8)
      ELTYPE='USED'
      GOTO 70                                                             20.01
  50  CALL MSGERR (3, 'Read error with variable '//NAAM_L)                40.03
      WRITE (PRINTF,41) NAAM, ELTEXT(1:LENCST)
      ELTYPE='USED'
!
  70  NS = LEN(CSTA)
      IF (NS.LT.LENW) THEN
        C(1:NS) = CSTA(1:NS)
        C(NS+1:LENW) = ' '
        LENCST = NS
      ELSE
        C(1:LENW) = CSTA(1:LENW)
        LENCST = LENW
      ENDIF
  80  IF (ITEST.GE.10) WRITE (PRINTF, 85) TRIM(NAAM), C, LENCST           40.03
  85  FORMAT (1X, A, ' = ', A, 4X, 'length:', I3)
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE INCTIM (IOPTIM, NAAM, RV, KONT, RSTA)
!                                                               *
!****************************************************************
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
!     30.72: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     30.72, Oct. 97: Introduced logical function EQREAL for floating point
!                     comparisons
!     30.04, Mar. 95: New subroutine
!     40.03, Feb. 00: local copy of NAAM used in error message
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Reads and interprets a time string
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     IOPTIM   int   inp   time reading option (see subr DTSTTI)
!
!
      INTEGER   IOPTIM
!
!     RV     : variable that is to be assigned a value
!     RSTA   : default value
!
      REAL      RV, RSTA
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
      CHARACTER NAAM *(*), KONT *(*)
!
!  5. PARAMETER VARIABLES
!
!     PARAMETERS: SEE PROGRAM DOCUMENTATION.
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     LENMN  : length of the string NAAM
!
      INTEGER    IENT, LENNM
!
!     NAAM_L : local copy of NAAM                                         40.03
!
      CHARACTER (LEN=40) :: NAAM_L                                        40.03
!
!     EQREAL : logical function, True if arguments are equal
!
      LOGICAL    EQREAL                                                   30.72
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE ( IENT, 'INCTIM')
!
      CHGVAL = .FALSE.
      NAAM_L = NAAM                                                       40.03
!     If necessary, a new data item is read.
      IF (ELTYPE.EQ.'USED') CALL  LEESEL
!     Consider type of data item.
      IF (ELTYPE.EQ.'KEY') THEN
!       find out whether in input is written:  NAAM=....
        LENNM = LEN(NAAM)
        IF (NAAM.NE.ELTEXT(1:LENNM)) GOTO 20
        ELTYPE = 'USED'
        CALL LEESEL
        IF (ELTYPE.EQ.'CHAR' .OR. ELTYPE.EQ.'OTHR'
     &    .OR. ELTYPE.EQ.'REAL' .OR. ELTYPE.EQ.'INT') GOTO 12
        GOTO 50
      ENDIF
  12  IF (ELTYPE.EQ.'CHAR' .OR. ELTYPE.EQ.'OTHR'
     &    .OR. ELTYPE.EQ.'REAL' .OR. ELTYPE.EQ.'INT') THEN
!       TYPE IS 'CHAR' etc., time is read from string
        CALL DTRETI (ELTEXT(1:LENCST), IOPTIM, RV)
        IF  (.NOT.EQREAL(RV,RSTA)) CHGVAL = .TRUE.                        30.72
        ELTYPE='USED'
        GOTO 80
      ENDIF
      IF (ELTYPE.EQ.'EOR') THEN
!       end of repeat group, standard value is assigned.
        IF (KONT.NE.'REP') GOTO 20
        ELTYPE='USED'
        GOTO 70
      ENDIF
      IF (ELTYPE.EQ.'EOF') GOTO 20
      IF (ELTYPE.EQ.'EMPT') THEN
        ELTYPE='USED'
        GOTO 20
      ENDIF
      IF (ELTYPE.EQ.'ERR') GOTO 50
      IF (ELTYPE.EQ.'INT') GOTO 40
      IF (ELTYPE.EQ.'REAL') GOTO 40
      CALL MSGERR (3, 'Error subr. INCTIM')                               40.00
      WRITE (PRINTF, '(1X,A,1X,A)') ELTYPE, KONT
      GOTO 80
!     --------------------------------------------------------
!     no value provided, action is taken according to 'KONT'.
!     --------------------------------------------------------
  20  IF (KONT.EQ.'REQ') GOTO 30
      IF (KONT.EQ.'REP') GOTO 70
      IF (KONT.EQ.'RQI') GOTO 28
      IF (KONT.EQ.'STA') GOTO 70
      IF (KONT.EQ.'NSKP') GOTO 70
      IF (KONT.EQ.'UNC') GOTO 80
      CALL MSGERR (3, 'Error subr. INCTIM')                               40.00
      WRITE (PRINTF, '(1X,A,1X,A)') ELTYPE, KONT
      GOTO 70                                                             20.01
!
  28  IF  (.NOT.EQREAL(RV,RSTA)) GOTO 80
  30  CALL MSGERR (3, 'No value for variable '//NAAM_L)                   40.03
      WRITE (PRINTF,31) NAAM, ELTEXT(1:LENCST)
  31  FORMAT (' -> ',A, '  item read=', A)
      GOTO 70                                                             20.01
  40  IF (KONT.EQ.'NSKP') GOTO 70
      CALL MSGERR (3, 'Wrong type of data for variable '//NAAM_L)         40.03
      WRITE (PRINTF,31) NAAM, ELTEXT(1:LENCST)
      ELTYPE='USED'
      GOTO 70                                                             20.01
  50  CALL MSGERR (3, 'Read error with variable '//NAAM_L)                40.03
      WRITE (PRINTF,31) NAAM, ELTEXT(1:LENCST)
      ELTYPE='USED'
!
  70  RV =RSTA
  80  IF (ITEST.GE.10) WRITE (PRINTF, 85) NAAM, ELTEXT(1:LENCST), RV
  85  FORMAT (1X, A, ' = ', A, 4X, 't in sec:', F10.0)
      RETURN
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE ININTV (NAME, RVAR, KONT, RSTA)                          30.09
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     Dec 1995, ver 30.09 : new subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Read a time interval in the form: number  DAY/HR/MIN/SEC
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     NAAM   : name of the variable according to the user manual
!     KONT   : What to do with the variable?
!              ='REQ'; error message if no value is found in the input file
!              ='UNC'; If no value, then variable will not be changed
!              ='STA'; If no value, then variable will get default value
!              ='RQI'; Variable may not have the value of RSTA
!              ='REP'  (repeat)
!              ='NSKP' (no skip) if data item is of different type,
!                      value is left unchanged.
!
      CHARACTER NAME *(*), KONT *(*)
!
!     RSTA   : default value
!     RVAR   : variable that is to be assigned a value
!
      REAL      RSTA, RVAR
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
!
!     FAC    : a factor, value depends on unit of time used
!     RI     : auxiliary variable
!
      REAL      FAC, RI
!
!     KEYWIS : logical function, True if keyword encountered is equal to
!              keyword in user manual
!
      LOGICAL   KEYWIS
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
!     -------------------------------------------------------------
!     Call INREAL to read number of time units
!     If a value was read
!     Then Read time unit
!          Case time unit is
!          DAY: Fac = 24*3600
!          HR:  Fac = 3600
!          MI:  Fac = 60
!          SEC: Fac = 1
!     Else Fac = 1
!     -------------------------------------------------------------
!     Interval in seconds = Fac * number of time units
!     -------------------------------------------------------------
!
! 13. SOURCE TEXT
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'ININTV')                                        30.09
!
      CALL INREAL (NAME, RI, KONT, RSTA)
      IF (CHGVAL) THEN
        CALL INKEYW ('STA', 'S')
        IF (KEYWIS('DA')) THEN
          FAC = 24.*3600.
        ELSE IF (KEYWIS('HR')) THEN
          FAC = 3600.
        ELSE IF (KEYWIS('MI')) THEN
          FAC = 60.
        ELSE
          CALL IGNORE ('S')
          FAC = 1.
        ENDIF
      ELSE
        FAC = 1.
      ENDIF
      RVAR = FAC * RI
      RETURN
!     end of subroutine ININTV                                            30.09
      END
!****************************************************************
!                                                               *
      SUBROUTINE LEESEL
!                                                               *
!****************************************************************
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
!     Jan. 1994, mod. 20.05: ELREAL is made double precision
!     40.13, Jan. 01: ! is now added as comment sign
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     reads a new data item from the string 'KAART'.
!     type of the item is determined, and the contents appears
!     in ELTEXT, ELINT, or ELREAL, as the case may be.
!     the following types are distinguished:
!     'KEY'   keyword
!     'INT'   integer or real number
!     'REAL'  real number
!     'CHAR'  character string enclosed in quotes
!     'EMPT'  empty data field
!     'OTHR'  non-empty data item not recognized as real, int or char,
!             possibly a time string
!     'EOF'   end of input file
!
!     'EOR'   end of repeat, or end of record
!     'ERR'   error
!     'USED'  used, item last read is processed already.
!
!  3. METHOD
!
!     difference between comment signs $ and !:                           40.13
!     everything on an input line behind a ! is ignored
!     text between two $-signs (on one line) is intepreted as comment
!     text behind two $-signs is intepreted as valid input
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     IRK    : auxiliary value used to detect errors
!     ISIGN1 : sign of mantissa part
!     ISIGN2 : sign of exponent part
!     ISTATE : state of the number reading process
!     J      : counter
!     JJ     : counter
!     JKAR   : counts the number of characters in the data field
!     NREP   : repetition number
!     NUM1   : value of integer part of mantissa
!     NUM2   : exponent value
!
      INTEGER   IENT, IRK, ISIGN1, ISIGN2, ISTATE, J, JJ, JKAR, NREP,
     &          NUM1, NUM2
!
!     RMANT  : real mantissa value
!
      DOUBLE    PRECISION RMANT
!
!     QUOTE  : the quote character
!
      CHARACTER QUOTE *1
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
      SAVE  IENT, QUOTE, NREP                                             40.00
      DATA  QUOTE/''''/ , IENT/0/, NREP/1/
      CALL  STRACE ( IENT, 'LEESEL')
!
      IF (NREP.GT.1) THEN
        NREP = NREP - 1
        GOTO 190
      ENDIF
!
!     initialisations
!
   2  NREP = 1
      DO  4 J=1,LINELN,4                                                  40.00
        ELTEXT(J:J+3) = '    '
   4  CONTINUE
      JKAR = 1
      ELINT=0
      ELREAL=0.
!
!     start processing data item
!
      IF (KARNR.EQ.0) GOTO 12
!     process a new character
  10  IF (KAR.EQ.'!' .OR. KARNR.GT.LINELN) THEN                           40.13
!       end of the line is reached, if repetition factor is >1
!       the data item is assumed to be empty
        IF (NREP.GT.1) GOTO 28
!       end of the line is reached, if no repetition factor appears
!       the data item is assumed to be of type 'EOR'
        ELTYPE='EOR'
        IF (KAR.EQ.'!') KARNR = LINELN+1                                  40.13
        GOTO 190
      ENDIF
!     skip leading blanks or Tab characters
  11  IF (KAR.NE.' ' .AND. KAR.NE.TABC) GOTO 20
  12  CALL GETKAR
!     end of input file was reached
      IF (ELTYPE.EQ.'EOF') THEN
!       generate keyword STOP
        ELTEXT='STOP'
        GOTO 190
      ENDIF
      GOTO 10
!     if character is comma, empty data field
  20  IF (KAR.NE.',') GOTO 30
      CALL GETKAR
  28  ELTYPE='EMPT'
      GOTO 190
!     Notice: jump to label 28 (empty data field)
!     if after repetition a comment, a keyword, end of record etc. is found.
!     --------------------------------------------------------
!     see whether end of repeat (; or /) is marked
  30  IF (INDEX(';/',KAR).GT.0) THEN                                      40.00
        IF (NREP.GT.1) GOTO 28
        ELTYPE='EOR'
        CALL GETKAR
        GOTO 190
      ENDIF
!     ( marks the beginning of a data item group; is ignored
  38  IF (KAR.EQ.'(') GOTO 12
!     --------------------------------------------------------
!     comment; data enclosed in comment identifiers is interpreted as comment
  40  IF (KAR.EQ.COMID) THEN
        IF (NREP.GT.1) GOTO 28
  41    CALL  GETKAR
        IF (KARNR.GT.LINELN) GOTO 10                                      40.00
        IF (KAR.NE.COMID) GOTO 41
        GOTO 12
      ENDIF
!     -------------------------------------------------------
!     if item is a number, read this integer or real number
!
!     integer number:  SIGN1]NUM1
!     real:            SIGN1]NUM1].]MANT]E]SIGN2]NUM2
!        ISTATE =    10     9    8 7    6 5     4    3
!     SIGN1, SIGN2:     + OR -
!     NUM1, NUM2, MANT: series of digits
!     -------------------------------------------------------
  50  IF (INDEX('+-.0123456789',KAR).EQ.0) GOTO 80
      NUM1=0
      NUM2=0
      ISIGN1=1
      ISIGN2=1
      ISTATE=10
      IRK=0
      RMANT=0.
      ELTYPE='INT'
      IF (INDEX('+-',KAR).EQ.0) GOTO 52
      ISTATE=9
      IF (KAR.EQ.'-') ISIGN1=-1
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
!     ****  part before decimal point  ****
  52  IF (INDEX('0123456789',KAR).EQ.0) GOTO 54
      IRK=1
      ISTATE=8
      NUM1=10*NUM1+INDEX('123456789',KAR)
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
      GOTO 52
  54  IF (KAR.NE.'.') GOTO 56
      ISTATE=7
      ELTYPE='REAL'
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
  56  JJ=-1
!     ****  part after decimal point  ****
  57  IF (INDEX('0123456789',KAR).EQ.0) GOTO 58
      IRK=1
      ISTATE=6
      RMANT = RMANT + DBLE(INDEX('123456789',KAR))*1.D1**JJ               20.05
      JJ=JJ-1
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
      GOTO 57
  58  IF (ISTATE.GE.9 .OR. IRK.EQ.0) GOTO 120
!     ****  exponent part  ****
      IF (INDEX('DdEe^',KAR).EQ.0) GOTO 66
      ISTATE=5
      IRK=0
      IF (ELTYPE.EQ.'INT') ELTYPE='REAL'
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
      IF (INDEX('+-',KAR).EQ.0) GOTO 62
      IF (KAR.EQ.'-') ISIGN2=-1
      ISTATE=4
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
  62  IF (INDEX('0123456789',KAR).EQ.0) GOTO 66
      IRK=1
      ISTATE=3
      NUM2=10*NUM2+INDEX('123456789',KAR)
      CALL PUTKAR (ELTEXT, KAR, JKAR)
      CALL GETKAR
      GOTO 62
!     ****  a number is put together  ****
  66  IF (IRK.EQ.0) GOTO 120
      IF (INDEX('+-.',KAR).GE.1) ELTYPE='OTHR'
      ISTATE=2
      IF (ITEST.GE.330) WRITE (PRINTF,699) ELTYPE, ISIGN1, NUM1,
     &  RMANT, ISIGN2, NUM2
 699  FORMAT (1X, A4, 2I6, F12.9, 2I6)
      IF (ELTYPE.EQ.'REAL') ELREAL =
     &  ISIGN1*(DBLE(NUM1)+RMANT) * 1.D1**(ISIGN2*NUM2)                   20.05
      IF (ELTYPE.EQ.'INT') ELINT = ISIGN1*NUM1
      LENCST = JKAR - 1                                                   30.03
!     skip trailing blanks
  67  IF (KAR.NE.' ' .AND. KAR.NE.TABC) GOTO 68
      ISTATE=1
      CALL GETKAR
      GOTO 67
!     If a * is encountered now, it is interpreted as a repetition factor.
  68  IF (KAR.EQ.'*') THEN
        IF (ELTYPE.EQ.'INT' .AND. ELINT.GT.0) THEN
          NREP = ELINT
          ELINT = 0
          CALL GETKAR
          GOTO 10
        ELSE
          CALL MSGERR (2, 'Wrong repetition factor')
          CALL GETKAR
          GOTO 190
        ENDIF
      ENDIF
  69  IF (KAR.EQ.',') THEN
        CALL GETKAR
        GOTO 190
      ENDIF
      IF (ISTATE.EQ.1) GOTO 190
      IF (INDEX(' ;',KAR).NE.0 .OR. KAR.EQ.TABC) THEN
        GOTO 190
      ENDIF
!     number is not followed by , blank or tab; type is made OTHR:
      GOTO 120
!     ----------------------------------------------------------
!     a character string is read; it start and ends with a quote
!     ----------------------------------------------------------
  80  IF (KAR.EQ.QUOTE) THEN                                              40.00
        ELTYPE='CHAR'
        LENCST = 0                                                          30.02
        JJ=1
  82    CALL GETKAR
!       end of the string: end of record or closing quote
        IF (KARNR.GT.LINELN) GOTO 190                                     40.00
        IF (KAR.EQ.QUOTE) THEN
          CALL GETKAR
!         new character is not a quote; end of the string
          IF (KAR.NE.QUOTE) GOTO 88
!         double quote is read as a single quote; continue
        ENDIF
!       put the character into ELTEXT
  84    ELTEXT(JJ:JJ) = KAR
        LENCST = JJ
        JJ=JJ+1
        GOTO 82
!       process characters behind the string
  87    CALL GETKAR
!       skip trailing blanks
  88    IF (KAR.EQ.' ' .OR. KAR.EQ.TABC) GOTO 87
        IF (KAR.NE.',') GOTO 190
        CALL GETKAR
        GOTO 190
      ENDIF
!     -------------------------------------------------------
!     a keyword is read
!     a keyword starts with a letter (upper or lower case)
!     -------------------------------------------------------
  90  CALL UPCASE (KAR)
      IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',KAR).GT.0) THEN              40.00
        IF (NREP.GT.1) GOTO 28
        ELTYPE='KEY'
        ISTATE=2
        JJ=1
  92    ELTEXT(JJ:JJ) = KAR
        LENCST = JJ                                                       30.02
        CALL GETKAR
        CALL UPCASE (KAR)
        JJ=JJ+1
!       next characters: letters, digits or - _ .
        IF (INDEX('ABCDEFGHIJKLMNOPQRSTUVWXYZ',KAR).GE.1) GOTO 92
        IF (INDEX('0123456789-_.',KAR).GE.1) GOTO 92
!       keyword is read
        KEYWRD = ELTEXT(1:8)
!       trailing blanks or tab char are skipped
  94    IF (KAR.NE.' ' .AND. KAR.NE.TABC) GOTO  96
        CALL GETKAR
        GOTO 94
!       closure character  : or = is processed
  96    IF (INDEX('=:',KAR).EQ.0) GOTO 190
        CALL GETKAR
        GOTO 190
      ENDIF
!     --------------------------------------------------
!     continuation symbol is read
!     --------------------------------------------------
 100  IF (INDEX('_&',KAR).EQ.0) GOTO 120
      IF (NREP.GT.1) GOTO 28
 110  KARNR=0
      GOTO 12
!     --------------------------------------------------
!     other type of data
!     --------------------------------------------------
 120  ELTYPE='OTHR'
 122  ELTEXT(JKAR:JKAR) = KAR                                             30.04
      LENCST = JKAR                                                       30.02
      JKAR=JKAR+1
      CALL GETKAR
      IF (INDEX(' ,;', KAR).GE.1 .OR. KAR.EQ.TABC) GOTO 126
      GOTO 122
 126  CALL GETKAR
! 127  CALL MSGERR (3, 'Read error in: ')
!      WRITE (PRINTF,129) ELTEXT
! 129  FORMAT (A)                                                         40.00
!      RETURN
!     --------------------------------------------------
!     test output and return to calling program
!     --------------------------------------------------
 190  IF (ITEST.GE.120) WRITE (PRTEST, 199) KAR, KARNR, ELTYPE, ELREAL,
     &  ELINT, NREP, ELTEXT(1:LENCST)
 199  FORMAT (' test LEESEL: ', A1, 1X, I4, 1X, A4, D12.4, 2I6, 2X, A)    20.05
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE GETKAR
!                                                               *
!****************************************************************
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
!     40.13, Jan. 2001: TRIM used to limit output
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This procedure reads a next character (KAR) from the string KAART.
!     The position of this character in KAART is indicated by KARNR.
!     If needed, a new input line is read.
!     At the end of the input file ELTYPE is made 'EOF'.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA  IENT /0/
      CALL STRACE (IENT, 'GETKAR')
      IF (KARNR.EQ.0) THEN
        READ (INPUTF, 7, END=20) KAART
   7    FORMAT (A)                                                        40.00
        IF (ITEST.GE.-10) WRITE (PRINTF, 8) TRIM(KAART)                   40.13
   8    FORMAT (1X,A)                                                     40.00
        KARNR=1
      ENDIF
      IF (KARNR.GT.LINELN) THEN                                           40.00
        KAR=';'
        GOTO 90
      ENDIF
      KAR = KAART(KARNR:KARNR)
      KARNR=KARNR+1
      GOTO 90
!     end of file is encountered
  20  ELTYPE='EOF'
      KAR='@'
  90  IF (ITEST.GE.320) WRITE (PRINTF, 99) ELTYPE, KAR, KARNR
  99  FORMAT (' Test GETKAR', 2X, A4, 2X, A1, I4)
      RETURN
!     end of subroutine GETKAR
      END
!****************************************************************
!                                                               *
      SUBROUTINE PUTKAR (LTEXT, KARR, JKAR)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     this procedure inserts a character (KARR) usually read by GETKAR
!     into the string LTEXT, usually equal to ELTEXT, in the place
!     JKAR. After this JKAR is increased by 1.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     JKAR   : counts the number of characters in a data field
!
      INTEGER  JKAR
!
!     LTEXT  : a character string; after a number of calls it should
!              contain the character representation of a data field
!     KARR   : character to be inserted into LTEXT
!
      CHARACTER LTEXT *(*), KARR *1
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'PUTKAR')
      IF (JKAR.GT.LEN(LTEXT)) CALL MSGERR (2, 'PUTKAR, string too long')
      LTEXT(JKAR:JKAR) = KARR
      LENCST = JKAR
      JKAR = JKAR + 1
      RETURN
!     end of subroutine PUTKAR
      END
!****************************************************************
!                                                               *
      SUBROUTINE UPCASE (CHARST)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     changes all characters of the string CHARST from lower to
!     upper case
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     CHARST : a character string
!
      CHARACTER*(*) CHARST                                                40.31
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IC     : sequence number of a character in the string CHARST
!     IENT   : Number of entries into this subroutine
!     KK     : position of a character in a given string
!     LLCC   : length of the given character string
!
      INTEGER   IC, IENT, KK, LLCC
!
!     ABCUP  : A to Z upper case characters
!     ABCLO  : a to z lower case characters
!     CC     : a character
!
      CHARACTER ABCUP *26, ABCLO *26, CC *1
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
      SAVE IENT
      DATA  IENT /0/
      DATA  ABCUP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA  ABCLO /'abcdefghijklmnopqrstuvwxyz'/
      CALL STRACE (IENT, 'UPCASE')
!
      LLCC = LEN (CHARST)
      DO 10 IC = 1, LLCC
         CC = CHARST(IC:IC)
         KK = INDEX (ABCLO, CC)
         IF (KK.NE.0) CHARST(IC:IC) = ABCUP(KK:KK)
  10  CONTINUE
      RETURN
!     end of subroutine UPCASE
      END
!****************************************************************
!                                                               *
      LOGICAL FUNCTION EQCSTR (STR1, STR2)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     EQCSTR is assigned the value true if the two strings are the same
!     (case-insensitive)
!
!  3. METHOD
!
!     both strings are converted to upper case and then compared
!     the length of string 2 gives the number of significant characters
!
!  4. ARGUMENT VARIABLES
!
      CHARACTER (LEN=*) :: STR1, STR2
!     two character strings to be compared
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
      INTEGER, SAVE  :: IENT = 0
!     IENT   : Number of entries into this subroutine

      INTEGER :: IC, LLCC
!     IC     : sequence number of a character in the string
!     LLCC   : length of the given character string

      CHARACTER (LEN=1) :: CC1, CC2
!     a character, one from STR1 the other from STR2
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
      CALL STRACE (IENT, 'UPCASE')
!
      EQCSTR = .TRUE.
      LLCC = LEN (STR2)
      IF (LEN(STR1).LT.LLCC) THEN
        EQCSTR = .FALSE.
        GOTO 90
      ENDIF
      DO IC = 1, LLCC
         CC1 = STR1(IC:IC)
         CALL UPCASE (CC1)
         CC2 = STR2(IC:IC)
         CALL UPCASE (CC2)
         IF (CC1.NE.CC2) THEN
           EQCSTR = .FALSE.
           GOTO 90
         ENDIF
      ENDDO
  90  RETURN
      END function EQCSTR
!****************************************************************
!                                                               *
      LOGICAL FUNCTION KEYWIS (STRING)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.00, July
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This procedure tests whether a keyword given by the user
!     coincides with a keyword known in the program (i.e. string).
!     if so, keywis is made .True., otherwise it is .False.
!     also ELTYPE is made 'USED', so that next element can be read.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     STRING : a keyword which is compared with a keyword found in the input file
!
      CHARACTER STRING *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!     J      : counter
!     LENSS  : length of the keyword STRING
!
      INTEGER   IENT, J, LENSS
!
!     KAR1   : a character of the keyword appearing in the input file
!     KAR2   : corresponding character in the STRING
!
      CHARACTER KAR1 *1, KAR2 *1
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'KEYWIS')
!
      KEYWIS = .FALSE.
      IF (ELTYPE.EQ.'USED') GOTO 30
!
      KEYWIS=.TRUE.
      LENSS = LEN (STRING)
      DO  20  J=1, LENSS
        KAR1 = KEYWRD(J:J)
        KAR2 = STRING(J:J)
        IF (KAR1.NE.KAR2 .AND. KAR2.NE.' ') THEN                          40.00
          KEYWIS=.FALSE.
          GOTO 30
        ENDIF
  20  CONTINUE
      IF (ELTYPE.EQ.'KEY') ELTYPE = 'USED'
  30  RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE  WRNKEY
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     THIS PROCEDURE PRODUCES AN ERROR MESSAGE
!     IT IS CALLED IF AN ILLEGAL KEYWORD IS FOUND IN THE
!     USER'S INPUT. IT MAKES ELTYPE = 'USED'
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'WRNKEY')
!
      CALL MSGERR (2, 'Illegal keyword: '//KEYWRD)
      ELTYPE = 'USED'
      RETURN
      END
!****************************************************************
!                                                               *
      SUBROUTINE  IGNORE (STRING)
!                                                               *
!****************************************************************
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
!     40.41: Marcel Zijlema
!
!  1. UPDATES
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     This procedure calls subroutine INKEYW to read a keyword.
!     if this keyword is equal to string, eltype is made 'USED'.
!     it is used if a keyword can occur in the input which
!     does not lead to any action.
!
!  3. METHOD
!
!  4. ARGUMENT VARIABLES
!
!     STRING : keyword (if appearing in input file) that can be ignored
!
      CHARACTER STRING *(*)
!
!  5. PARAMETER VARIABLES
!
!  6. LOCAL VARIABLES
!
!     IENT   : Number of entries into this subroutine
!
      INTEGER   IENT
!
!     KEYWIS : logical function
!
      LOGICAL   KEYWIS
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
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'IGNORE')
!
      CALL INKEYW ('STA', 'XXXX')
      IF (KEYWIS(STRING)) RETURN
      IF (KEYWIS('XXXX')) RETURN
      IF (ITEST.GE.60) WRITE (PRINTF, 5) KEYWRD, ELTYPE
   5  FORMAT (' NOT IGNORED: ', A, 2X, A)
      RETURN
      END
