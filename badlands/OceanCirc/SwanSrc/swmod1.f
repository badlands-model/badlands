!                 COMMON VARIABLES RELATED MODULES, file 1 of 3
!
!     Contents of this file
!
!     OCPCOMM1           contains common variables for Ocean Pack
!     OCPCOMM2           contains common variables for Ocean Pack
!     OCPCOMM3           contains common variables for Ocean Pack
!     OCPCOMM4           contains common variables for Ocean Pack
!     SWCOMM1            contains common variables for SWAN
!     SWCOMM2            contains common variables for SWAN
!     SWCOMM3            contains common variables for SWAN
!     SWCOMM4            contains common variables for SWAN
!     TIMECOMM           contains common variables for SWAN
!
      MODULE OCPCOMM1
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM1.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     LINELN [ 120]  max length of input lines (in command file)
!
      INTEGER LINELN
      PARAMETER  (LINELN=120)
!
!  7. Local variables
!
!     *** character data used by the command reading system ***
!
! BLANK  ['    '] blank string
! COMID  [   '$'] character which distinguishes comments in the command input
! ELTEXT [      ] contents of the last string read by reading system
! ELTYPE [      ] type of the element last read by reading system
!                 ='CHAR'; last read data element is the string in ELTEXT
!                 ='EMPT;  empty data field
!                 ='EOF';  end of file has been reached
!                 ='EOR';  end of repeat has been reached
!                 ='ERR';  incorrect data field was encountered
!                 ='INT';  last read data element is the integer in ELINT
!                 ='KEY';  last read data element is the keyword in KEYWRD
!                 ='OTHR'; other
!                 ='REAL'; last read data element is the real in ELREAL
!                 ='USED'; last read data element is processed, new can be read
! KAART  [      ] contents of the input line last read by reading system
! KAR    [      ] character last read by reading system
!                 =COMID; begin or end of comment
!                 =TABC; data element separation mark
!                 =' '; data element separation mark
!                 ='&'; continuation mark
!                 ='('; begin of data group
!                 ='*'; repetition mark
!                 =','; data element seperation mark (only for numbers)
!                 ='/'; end of repetition mark
!                 =':'; assignment mark
!                 =';'; end of record or input line
!                 ='='; assignment mark
!                 ='@'; end of file mark
!                 ='_'; continuation mark
!                 other: letter or digit to be processed by reading system
! KEYWRD [      ] contents of the last keyword read by reading system
! TABC   [CALCUL] =CHAR(9); tabular character
!
      CHARACTER*4        BLANK
      CHARACTER          COMID
      CHARACTER*(LINELN) ELTEXT
      CHARACTER*4        ELTYPE
      CHARACTER*(LINELN) KAART
      CHARACTER          KAR
      CHARACTER*8        KEYWRD
      CHARACTER          TABC
!
!     *** numerical data used by the command reading system ***
!
! CHGVAL [      ] whether last read value is different from a given value for
!                 subroutines INREAL, ININTG, INCSTR, INCTIM
! ELINT  [      ] last element read from user command, when integer
! KARNR  [      ] position on the input line of character last processed
!                 by the reading system,
!                 =0; no characters read yet
!                 =81; next input line has to be read to the common KAART first
! LENCST [      ] length of the string stored in ELTEXT
! ELREAL [      ] last element read from user command, when real or double
!
      INTEGER          ELINT, KARNR, LENCST
      DOUBLE PRECISION ELREAL
      LOGICAL          CHGVAL
!
!     *** origin for day and time ***
!
! REFDAY    [    ] Day number of the reference day. The first day entered is used
!                  as reference day, the reference time is 0:00 of the reference day.
!
      INTEGER REFDAY
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM1

      MODULE OCPCOMM2
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM2.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     LENFNM [  80]  length of file names
!
      INTEGER LENFNM
      PARAMETER (LENFNM=80)
!
!  7. Local variables
!
!     *** names and other character strings ***
!
! DIRCH1 [\     ] directory separation character as appears in input file
! DIRCH2 [\     ] directory separation character replacing DIRCH1
! FILEA  [      ] not used
! FILEB  [      ] not used
! FILENM [      ] file name of the file currently used for I/O
! INST   ['Delft University of Technology'] name of the institute
!                 Can be changed in the file SWANINIT
! PROJID ['SWAN'] acronym of the project for which the computation is taking place
!                 ='NAME'; set by command PROJ 'NAME' ...
! PROJNR [CALCUL] =BLANK; run number for the computation
!                 ='NR'; set by command PROJ ... 'NR' ...
! PROJT1 [CALCUL] =BLANK; 1st line of the project title
!                 ='title1'; set by command PROJ ... 'title1' ...
! PROJT2 [CALCUL] =BLANK; 2nd line of the project title
!                 ='title2'; set by command PROJ ... 'title2' ...
! PROJT3 [CALCUL] =BLANK; 3rd line of the project title
!                 ='title3'; set by command PROJ ... 'title3'
! PTITLE [      ] not used
! VERTXT [calcul] program version, character representation
!
      CHARACTER (LEN=1)      :: DIRCH1, DIRCH2
      CHARACTER (LEN=LENFNM) :: FILEA , FILEB , FILENM
      CHARACTER (LEN=40)     :: INST
      CHARACTER (LEN=16)     :: PROJID
      CHARACTER (LEN=4)      :: PROJNR
      CHARACTER (LEN=72)     :: PROJT1, PROJT2, PROJT3
      CHARACTER (LEN=36)     :: PTITLE
      CHARACTER (LEN=20)     :: VERTXT
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM2

      MODULE OCPCOMM3
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM3.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     ---
!
!  7. Local variables
!
!     *** data for output, mainly plotting ***
!
! DXQ    [      ] mesh size of the output frame in X-direction
!                 =0.01; if MXQ=1
!                 =XQLEN/(MXQ-1); if MXQ>1
! DYQ    [      ] mesh size of the output frame in Y-direction
!                 =0.01; if MYQ=1
!                 =YQLEN/(MYQ-1); if MYQ>1
! MXQ    [CALCUL] number of grid points of the output frame in X-direction
! MYQ    [CALCUL] number of grid points of the output frame in Y-direction
! VERNUM [ 40.41] version number of SWAN
!
      INTEGER MXQ, MYQ
      REAL    DXQ, DYQ, VERNUM
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM3

      MODULE OCPCOMM4
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM4.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     INAN   [   1]  integer representing not a number
!     RNAN   [   1]  real representing not a number
!
      INTEGER INAN
      REAL    RNAN
      PARAMETER (INAN=-1073750760,
     &           RNAN=-1.07374515E+09)
!
!  7. Local variables
!
!     *** file unit reference numbers ***
!
! EXPORT [   ] not used
! FUNHI  [CAL] highest free unit number
!              =IUNMAX
! FUNLO  [ 21] lowest free unit number
! HIOPEN [   ] highest unit number of an open file
! IMPORT [   ] not used
! INPUTF [  3] unit number for the file with command input ('INPUT')
!              set by file SWANINIT
! ITMOPT [  1] time coding option
!              =1; ????
!              =3; ????
!              set by file SWANINIT
! IUNMIN [  0] minimum unit number
! IUNMAX [ 99] maximum unit number
!              set by file SWANINIT
! PRINTF [  4] unit number for the file with standard output ('PRINT')
!              set by file SWANINIT
! PRTEST [CAL] unit number for the print file containing test output
!              =PRINTF
!              set by file SWANINIT
! SCREEN [  6] unit number for the screen
!              (is for batch-oriented systems equal to PRINTF)
!
      INTEGER EXPORT, FUNHI , FUNLO , HIOPEN
      INTEGER IMPORT, INPUTF, ITMOPT, IUNMAX
      INTEGER IUNMIN, PRINTF, PRTEST, SCREEN

!     *** test parameters ***
!
! ITEST  [  0] indicates the amount of test output requested
!              =30; for command TEST
!              =itest; set by command TEST [itest] [itrace]
! ITRACE [  0] a message is printed up to ITRACE times
!              =itrace; set by command TEST [itest] [itrace]
! LEVERR [  0] severity of the errors encountered.
! LTRACE [.F.] indicates whether to call STRACE
!              =.T.; when ITRACE>0
! MAXERR [  1] maximum severity of errors allowed, if larger no computation
!              =1; warnings
!              =2; errors
!              =3; severe errors
!              =4; terminating errors
!              =maxerr; set by command SET ... [maxerr] ...
!
      INTEGER ITEST, ITRACE, LEVERR, MAXERR
      LOGICAL LTRACE
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM4

      MODULE SWCOMM1
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
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!     40.61: Marcel Zijlema
!     40.64: Marcel Zijlema
!     41.12: Nico Booij
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file SWCOMM1.INC
!     40.21, Apr. 04: NMOVAR increased by 1
!                     (reserved for diffraction)
!     40.41, Sep. 04: NMOVAR increased by 2 (quantities TMM10
!                     and RTMM10 added)
!     40.51, Feb. 05: NMOVAR increased by 1 (quantity TMBOT added)
!     40.51, Sep. 05: NMOVAR increased by 2 (quantities WATLEV
!                     and BOTLEV added)
!     40.51, Feb. 06: NMOVAR increased by 1 (quantity TPS added)
!     40.61, Sep. 06: NMOVAR increased by 3 (quantities DISBOT,
!                     DISSRF and DISWCP added)
!     40.61, Sep. 06: NMOVAR increased by 1 (quantity DISVEG added)
!     40.64, Apr. 07: NMOVAR increased by 2 (quantities QP and BFI added)
!     40.85, Aug. 08: NMOVAR increased by 10 (quantities GENE, GENW, REDI, REDQ, REDT,
!                                             PROPA, PROPX, PROPT, PROPS and RADS added)
!     41.12, Apr. 10: NMOVAR increased by 1 (quantity NPL added)
!     41.15, Mar. 11: NMOVAR increased by 1 (quantity LWAVP added)
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     NMOVAR [  70]  maximum number of output variables                   40.85 40.64 40.61 40.51 40.21 40.41
!     MOUTPA [  50]  number of output parameters                          40.87
!
      INTEGER NMOVAR, MOUTPA
      PARAMETER (NMOVAR = 71, MOUTPA=50)                                  41.15 40.85 40.87 40.64 40.61 40.51 40.21 40.41
!
!  7. Local variables
!
!     *** names and other character data ***
!
! CHTIME [ '    '] character string representation of date-time of computation
! FNEST [CALCULAT] BCF (=BLANK), name of nest file
! OVKEYW(NMOVAR)   keyword identifying output quantity in a SWAN command
!  ( 1) [    'XP']
!  ( 2) [    'YP']
!  ( 3) [  'DIST']
!  ( 4) [   'DEP']
!  ( 5) [   'VEL']
!  ( 6) [  'UBOT']
!  ( 7) [  'DISS']
!  ( 8) [    'QB']
!  ( 9) [   'LEA']
!  (10) [    'HS']
!  (11) [  'TM01']
!  (12) [   'RTP']
!  (13) [   'DIR']
!  (14) [   'PDI']
!  (15) [   'TDI']
!  (16) [  'DSPR']
!  (17) [  'WLEN']
!  (18) [  'STEE']
!  (19) [   'TRA']
!  (20) [   'FOR']
!  (21) ['      ']
!  (22) ['      ']
!  (23) ['      ']
!  (24) [    'XC']
!  (25) [    'YC']
!  (26) [  'WIND']
!  (27) [   'FRC']
!  (28) [ 'RTM01']
!  (29) ['      ']
!  (30) [   'DHS']
!  (31) ['DRTM01']
!  (32) [  'TM02']
!  (33) [  'FSPR']
!  (34) [  'URMS']
!  (35) [  'UFRI']
!  (36) [  'ZLEN']
!  (37) [  'TAUW']
!  (38) [ 'CDRAG']
!  (39) [ 'SETUP']
!  (40) [  'TIME']
!  (41) [  'TSEC']
!  (47) [ 'TMM10']
!  (48) ['RTMM10']
!  (49) ['DIFPAR')
!  (50) [ 'TMBOT']
!  (51) [  'WATL']
!  (52) [  'BOTL']
!  (53) [   'TPS']
!  (54) [  'DISB']
!  (55) [ 'DISSU']
!  (56) [  'DISW']
!  (57) [  'DISV']
!  (58) [    'QP']
!  (59) [   'BFI']
!  (60) [  'GENE']
!  (61) [  'GENW']
!  (62) [  'REDI']
!  (63) [  'REDQ']
!  (64) [  'REDT']
!  (65) [ 'PROPA']
!  (66) [ 'PROPX']
!  (67) [ 'PROPT']
!  (68) [ 'PROPS']
!  (69) [  'RADS']
!  (70) [   'NPL']
!  (71) [ 'LWAVP']
! OVSNAM(NMOVAR)   short name of output quantity
!  ( 1) [    'Xp']
!  ( 2) [    'Yp']
!  ( 3) [  'Dist']
!  ( 4) [ 'Depth']
!  ( 5) [   'Vel']
!  ( 6) [  'Ubot']
!  ( 7) ['Dissip']
!  ( 8) [    'Qb']
!  ( 9) [  'Leak']
!  (10) [    'Hs']
!  (11) [  'Tm01']
!  (12) [ 'Tpeak']
!  (13) [   'Dir']
!  (14) [ 'PkDir']
!  (15) [  'TDir']
!  (16) [  'Dspr']
!  (17) [  'Wlen']
!  (18) ['Steepn']
!  (19) ['Transp']
!  (20) ['WForce']
!  (21) ['AcDens']
!  (22) ['EnDens']
!  (23) [   'Aux']
!  (24) [    'Xc']
!  (25) [    'Yc']
!  (26) [ 'Windv']
!  (27) ['FrCoef']
!  (28) [ 'RTm01']
!  (29) ['EnDens']
!  (30) [   'dHs']
!  (31) [   'dTm']
!  (32) [  'Tm02']
!  (33) [  'FSpr']
!  (34) [  'Urms']
!  (35) [ 'Ufric']
!  (36) [  'Zlen']
!  (37) [  'TauW']
!  (38) [ 'Cdrag']
!  (39) [ 'Setup']
!  (40) [  'Time']
!  (41) [  'Tsec']
!  (47) [ 'Tm_10']
!  (48) ['RTm_10']
!  (49) ['DifPar']
!  (50) [ 'TmBot']
!  (51) ['Watlev']
!  (52) ['Botlev']
!  (53) ['TPsmoo']
!  (54) [ 'Sfric']
!  (55) [ 'Ssurf']
!  (56) [ 'Swcap']
!  (57) [  'Sveg']
!  (58) [    'Qp']
!  (59) [   'BFI']
!  (60) ['Genera']
!  (61) [ 'Swind']
!  (62) ['Redist']
!  (63) [  'Snl4']
!  (64) [  'Snl3']
!  (65) ['Propag']
!  (66) ['Propxy']
!  (67) ['Propth']
!  (68) ['Propsi']
!  (69) ['Radstr']
!  (70) ['Nplant']
!  (71) [ 'Lwavp']
! OVLNAM(NMOVAR)   long name of output quantity
!  ( 1) [                        'X user coordinate']
!  ( 2) [                        'Y user coordinate']
!  ( 3) [              'distance along output curve']
!  ( 4) [                                    'Depth']
!  ( 5) [                         'Current velocity']
!  ( 6) [           'Orbital velocity at the bottom']
!  ( 7) [                       'Energy dissipation']
!  ( 8) [                  'Fraction breaking waves']
!  ( 9) [     'Energy leak over spectral boundaries']
!  (10) [                  'Significant wave height']
!  (11) [             'Average absolute wave period']
!  (12) [                              'Peak period']
!  (13) [                   'Average wave direction']
!  (14) [    'direction of the peak of the spectrum']
!  (15) [        'direction of the energy transport']
!  (16) [                    'directional spreading']
!  (17) [                      'Average wave length']
!  (18) [                           'Wave steepness']
!  (19) [                    'Wave energy transport']
!  (20) [       'Wave driven force per unit surface']
!  (21) [                  'spectral action density']
!  (22) [                  'spectral energy density']
!  (23) [                       'auxiliary variable']
!  (24) [          'X computational grid coordinate']
!  (25) [          'Y computational grid coordinate']
!  (26) [    'Wind velocity at 10 m above sea level']
!  (27) [              'Bottom friction coefficient']
!  (28) [             'Average relative wave period']
!  (29) [ 'energy density integrated over direction']
!  (30) [      'difference in Hs between iterations']
!  (31) [      'difference in Tm between iterations']
!  (32) [                     'Zero-crossing period']
!  (33) [         'Frequency spectral width (Kappa)']
!  (34) [    'RMS of orbital velocity at the bottom']
!  (35) [                        'Friction velocity']
!  (40) [                                'Date-time']
!  (41) [      'Time in seconds from reference time']
!  (36) ['Zero velocity thickness of boundary layer']
!  (37) [                                     '    ']
!  (38) [                         'Drag coefficient']
!  (39) [                       'Setup due to waves']
!  (47) [             'Average absolute wave period']
!  (48) [             'Average relative wave period']
!  (49) [                    'Diffraction parameter']
!  (50) [                       'Bottom wave period']
!  (51) [                              'Water level']
!  (52) [                             'Bottom level']
!  (53) [            'Relative peak period (smooth)']
!  (54) [              'Bottom friction dissipation']
!  (55) [                'Surf breaking dissipation']
!  (56) [                 'Whitecapping dissipation']
!  (57) [                   'Vegetation dissipation']
!  (58) [                               'Peakedness']
!  (59) [                      'Benjamin-Feir index']
!  (60) [                        'Energy generation']
!  (61) [                         'Wind source term']
!  (62) [                    'Energy redistribution']
!  (63) [        'Total absolute 4-wave interaction']
!  (64) [        'Total absolute 3-wave interaction']
!  (65) [                       'Energy propagation']
!  (66) [                           'xy-propagation']
!  (67) [                        'theta-propagation']
!  (68) [                        'sigma-propagation']
!  (69) [                         'Radiation stress']
!  (70) [                            'Plants per m2']
!  (71) [                         'Peak wave length']
! OVUNIT(NMOVAR)   unit of of output quantity
!  ( 1) [CALCULAT] =UL
!  ( 2) [CALCULAT] =UL
!  ( 3) [CALCULAT] =UL
!  ( 4) [CALCULAT] =UH
!  ( 5) [CALCULAT] =UV
!  ( 6) [CALCULAT] =UV
!  ( 7) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  ( 8) [     ' ']
!  ( 9) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (10) [CALCULAT] =UH
!  (11) [CALCULAT] =UT
!  (12) [CALCULAT] =UT
!  (13) [CALCULAT] =UDI
!  (14) [CALCULAT] =UDI
!  (15) [CALCULAT] =UDI
!  (16) [CALCULAT] =UDI
!  (17) [CALCULAT] =UL
!  (18) [     ' ']
!  (19) [   'm2s'] changed to 'W/m' (=UP), if INRHOG=1
!  (20) [CALCULAT] =UF
!  (21) [   'm2s'] changed to 'Js/m2', if INRHOG=1
!  (22) [    'm2'] changed to 'J/m2', if INRHOG=1
!  (23) [     ' ']
!  (24) [     ' ']
!  (25) [     ' ']
!  (26) [CALCULAT] =UV
!  (27) [     ' ']
!  (28) [CALCULAT] =UT
!  (29) [    'm2'] changed to 'J/m2', if INRHOG=1
!  (30) [CALCULAT] =UH
!  (31) [CALCULAT] =UT
!  (32) [CALCULAT] =UT
!  (33) [     ' ']
!  (34) [CALCULAT] =UV
!  (35) [CALCULAT] =UV
!  (36) [CALCULAT] =UL
!  (37) [     ' ']
!  (38) [     ' ']
!  (39) [     'm']
!  (40) [     ' ']
!  (41) [     's']
!  (47) [CALCULAT] =UT
!  (48) [CALCULAT] =UT
!  (49) [     ' ']
!  (50) [CALCULAT] =UT
!  (51) [CALCULAT] =UH
!  (52) [CALCULAT] =UH
!  (53) [CALCULAT] =UT
!  (54) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (55) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (56) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (57) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (58) [     ' ']
!  (59) [     ' ']
!  (60) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (61) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (62) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (63) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (64) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (65) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (66) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (67) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (68) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (69) [  'm2/s'] change to UDL ???, changed to 'W/m2', if INRHOG=1
!  (70) [  '1/m2']
!  (71) [CALCULAT] =UL
! SNAME [        ] name of output point set
! UAP   [  'W/m2'] unit of dissipation
! UD    [        ] not used
! UDI   [  'degr'] unit of direction
! UDL   [  'm2/s'] unit of dissipation
! UET   [  'm3/s'] unit of energy transport, and wave force
! UF    [  'N/m2'] unit of pressure or shear stress (force per area)
! UH    [     'm'] unit of vertical length
! UL    [     'm'] unit of horizontal length
! UP    [   'W/m'] unit of energy flux density
! UST   [ 'm2/s2'] not used
! UT    [   'sec'] unit of time (change to s ???)
! UV    [   'm/s'] unit of velocity
!
      CHARACTER*20 CHTIME
      CHARACTER*36 FBCL,          FBCR,        FNEST
      CHARACTER*8  OVKEYW(NMOVAR)
      CHARACTER*40 OVLNAM(NMOVAR)
      CHARACTER*6  OVSNAM(NMOVAR)
      CHARACTER*16 OVUNIT(NMOVAR)
      CHARACTER*8  SNAME
      CHARACTER*6  UAP,           UD,          UDI,         UDL
      CHARACTER*6  UET,           UF,          UH,          UL
      CHARACTER*6  UP,            UST,         UT,          UV

!     *** information for output ***
!
! AKPOWR [       ] power in expression for computation of average wave number
!                  =kpower; set by command SET ... [kpower] ... (not documented)
! ALCQ   [       ] angle between x-axes of computational grid and output frame
! ALPQ   [       ] angle between x-axes of user coord. system and output frame
! COSCQ  [       ] cos of ALCQ
! COSPQ  [       ] cos of ALPQ
! DXK    [       ] mesh size of output frame
! DYK    [       ] mesh size of output frame
! ERRPTS [       ] unit ref. number of file containing coord. of "problem points"
!                  =16, unit reference number of the file
! INRHOG [      0] indicates the choice for output based on "variance" or "true energy"
!                  =0, output based on variance
!                  =1, output based on true energy
! IUBOTR [      0] set to 1, when IVTYPE=6 or 18
! OVSVTY(NMOVAR) type of the output variable
!                =1, scalar
!                =2, angle
!                =3, vector
!                =4, tensor
!                =5, fully spectral quantity
!                =6, directional spectral quantity
!  ( 1)  [  1]   Xp
!  ( 2)  [  1]   Yp
!  ( 3)  [  1]   Dist
!  ( 4)  [  1]   Depth
!  ( 5)  [  3]   Vel
!  ( 6)  [  1]   Ubot
!  ( 7)  [  1]   Dissip
!  ( 8)  [  1]   Qb
!  ( 9)  [  1]   Leak
!  (10)  [  1]   Hs
!  (11)  [  1]   Tm01
!  (12)  [  1]   Tpeak
!  (13)  [  2]   Dir
!  (14)  [  2]   PkDir
!  (15)  [  2]   TDir
!  (16)  [  1]   Dspr
!  (17)  [  1]   Wlen
!  (18)  [  1]   Steepn
!  (19)  [  3]   Transp
!  (20)  [  3]   WForce
!  (21)  [  5]   AcDens
!  (22)  [  5]   EnDens
!  (23)  [  1]   Aux
!  (24)  [  1]   Xc
!  (25)  [  1]   Yc
!  (26)  [  3]   Windv
!  (27)  [  1]   FrCoef
!  (28)  [  1]   RTm01
!  (29)  [  5]   EnDens
!  (30)  [  1]   dHs
!  (31)  [  1]   dTm
!  (32)  [  1]   Tm02
!  (33)  [  1]   FSpr
!  (34)  [  1]   Urms
!  (35)  [  1]   Ufric
!  (36)  [  1]   Zlen
!  (37)  [  1]   TauW
!  (38)  [  1]   Cdrag
!  (39)  [  1]   Setup
!  (47)  [  1]   Tm_10
!  (48)  [  1]   RTm_10
!  (49)  [  1]   DifPar
!  (50)  [  1]   TmBot
!  (51)  [  1]   Watlev
!  (52)  [  1]   Botlev
!  (53)  [  1]   TPSmoo
!  (54)  [  1]   Sfric
!  (55)  [  1]   Ssurf
!  (56)  [  1]   Swcap
!  (57)  [  1]   Sveg
!  (58)  [  1]   Qp
!  (59)  [  1]   BFI
!  (60)  [  1]   Genera
!  (61)  [  1]   Swind
!  (62)  [  1]   Redist
!  (63)  [  1]   Snl4
!  (64)  [  1]   Snl3
!  (65)  [  1]   Propag
!  (66)  [  1]   Propxy
!  (67)  [  1]   Propth
!  (68)  [  1]   Propsi
!  (69)  [  1]   Radstr
!  (70)  [  1]   Nplant
!  (71)  [  1]   Lwavp
! OVEXCV(NMOVAR)   exception value for output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [  -999.] Dist
!  ( 4)  [   -99.] Depth
!  ( 5)  [     0.] Vel
!  ( 6)  [   -10.] Ubot
!  ( 7)  [   -10.] Dissip
!  ( 8)  [    -1.] Qb
!  ( 9)  [   -10.] Leak
!  (10)  [   -10.] Hs
!  (11)  [   -10.] Tm01
!  (12)  [   -10.] Tpeak
!  (13)  [  -999.] Dir
!  (14)  [  -999.] PkDir
!  (15)  [  -999.] TDir
!  (16)  [   -10.] Dspr
!  (17)  [   -10.] Wlen
!  (18)  [    -1.] Steepn
!  (19)  [     0.] Transp
!  (20)  [     0.] WForce
!  (21)  [  -100.] AcDens
!  (22)  [  -100.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [ -1.E10] Xc
!  (25)  [ -1.E10] Yc
!  (26)  [     0.] Windv
!  (27)  [    -1.] FrCoef
!  (28)  [   -10.] RTm01
!  (29)  [  -100.] EnDens
!  (30)  [   -10.] dHs
!  (31)  [   -10.] dTm
!  (32)  [   -10.] Tm02
!  (33)  [    -1.] FSpr
!  (34)  [   -10.] Urms
!  (35)  [   -10.] Ufric
!  (36)  [    -1.] Zlen
!  (37)  [   -10.] TauW
!  (38)  [    -1.] Cdrag
!  (39)  [    -9.] Setup
!  (47)  [    -9.] Tm_10
!  (48)  [    -9.] RTm_10
!  (49)  [   -99.] DifPar
!  (50)  [    -9.] TmBot
!  (51)  [   -99.] Watlev
!  (52)  [   -99.] Botlev
!  (53)  [   -10.] TPSmoo
!  (54)  [    -9.] Sfric
!  (55)  [    -9.] Ssurf
!  (56)  [    -9.] Swcap
!  (57)  [    -9.] Sveg
!  (58)  [    -9.] Qp
!  (59)  [    -9.] BFI
!  (60)  [    -9.] Genera
!  (61)  [    -9.] Swind
!  (62)  [    -9.] Redist
!  (63)  [    -9.] Snl4
!  (64)  [    -9.] Snl3
!  (65)  [    -9.] Propag
!  (66)  [    -9.] Propxy
!  (67)  [    -9.] Propth
!  (68)  [    -9.] Propsi
!  (69)  [    -9.] Radstr
!  (70)  [    -9.] Nplant
!  (71)  [    -9.] Lwavp
! OVLEXP(NMOVAR)   lower expected limit of output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [     0.] Dist
!  ( 4)  [  -100.] Depth
!  ( 5)  [    -2.] Vel
!  ( 6)  [     0.] Ubot
!  ( 7)  [     0.] Dissip
!  ( 8)  [     0.] Qb
!  ( 9)  [     0.] Leak
!  (10)  [     0.] Hs
!  (11)  [     0.] Tm01
!  (12)  [     0.] Tpeak
!  (13)  [     0.] Dir
!  (14)  [     0.] PkDir
!  (15)  [     0.] TDir
!  (16)  [     0.] Dspr
!  (17)  [     0.] Wlen
!  (18)  [     0.] Steepn
!  (19)  [   -10.] Transp
!  (20)  [   -10.] WForce
!  (21)  [     0.] AcDens
!  (22)  [     0.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [  -100.] Xc
!  (25)  [  -100.] Yc
!  (26)  [   -50.] Windv
!  (27)  [     0.] FrCoef
!  (28)  [     0.] RTm01
!  (29)  [     0.] EnDens
!  (30)  [     0.] dHs
!  (31)  [     0.] dTm
!  (32)  [     0.] Tm02
!  (33)  [     0.] FSpr
!  (34)  [     0.] Urms
!  (35)  [     0.] Ufric
!  (36)  [     0.] Zlen
!  (37)  [     0.] TauW
!  (38)  [     0.] Cdrag
!  (39)  [    -1.] Setup
!  (47)  [     0.] Tm_10
!  (48)  [     0.] RTm_10
!  (49)  [   -10.] DifPar
!  (50)  [     0.] TmBot
!  (51)  [  -100.] Watlev
!  (52)  [  -100.] Botlev
!  (53)  [     0.] TPSmoo
!  (54)  [     0.] Sfric
!  (55)  [     0.] Ssurf
!  (56)  [     0.] Swcap
!  (57)  [     0.] Sveg
!  (58)  [     0.] Qp
!  (59)  [     0.] BFI
!  (60)  [     0.] Genera
!  (61)  [     0.] Swind
!  (62)  [     0.] Redist
!  (63)  [     0.] Snl4
!  (64)  [     0.] Snl3
!  (65)  [     0.] Propag
!  (66)  [     0.] Propxy
!  (67)  [     0.] Propth
!  (68)  [     0.] Propsi
!  (69)  [     0.] Radstr
!  (70)  [     0.] Nplant
!  (71)  [     0.] Lwavp
! OVLLIM(NMOVAR)   lower limit of validity of output quantity
!  ( 1)  [ -1.E10] Xp
!  ( 2)  [ -1.E10] Yp
!  ( 3)  [     0.] Dist
!  ( 4)  [  -1.E4] Depth
!  ( 5)  [  -100.] Vel
!  ( 6)  [     0.] Ubot
!  ( 7)  [     0.] Dissip
!  ( 8)  [     0.] Qb
!  ( 9)  [     0.] Leak
!  (10)  [     0.] Hs
!  (11)  [     0.] Tm01
!  (12)  [     0.] Tpeak
!  (13)  [     0.] Dir
!  (14)  [     0.] PkDir
!  (15)  [     0.] TDir
!  (16)  [     0.] Dspr
!  (17)  [     0.] Wlen
!  (18)  [     0.] Steepn
!  (19)  [  -100.] Transp
!  (20)  [  -1.E5] WForce
!  (21)  [     0.] AcDens
!  (22)  [     0.] EnDens
!  (23)  [ -1.E10] Aux
!  (24)  [ -1000.] Xc
!  (25)  [ -1000.] Yc
!  (26)  [  -100.] Windv
!  (27)  [     0.] FrCoef
!  (28)  [     0.] RTm01
!  (29)  [     0.] EnDens
!  (30)  [     0.] dHs
!  (31)  [     0.] dTm
!  (32)  [     0.] Tm02
!  (33)  [     0.] FSpr
!  (34)  [     0.] Urms
!  (35)  [     0.] Ufric
!  (36)  [     0.] Zlen
!  (37)  [     0.] TauW
!  (38)  [     0.] Cdrag
!  (39)  [    -1.] Setup
!  (47)  [     0.] Tm_10
!  (48)  [     0.] RTm_10
!  (49)  [   -50.] DifPar
!  (50)  [     0.] TmBot
!  (51)  [  -1.E4] Watlev
!  (52)  [  -1.E4] Botlev
!  (53)  [     0.] TPSmoo
!  (54)  [     0.] Sfric
!  (55)  [     0.] Ssurf
!  (56)  [     0.] Swcap
!  (57)  [     0.] Sveg
!  (58)  [     0.] Qp
!  (59)  [     0.] BFI
!  (60)  [     0.] Genera
!  (61)  [     0.] Swind
!  (62)  [     0.] Redist
!  (63)  [     0.] Snl4
!  (64)  [     0.] Snl3
!  (65)  [     0.] Propag
!  (66)  [     0.] Propxy
!  (67)  [     0.] Propth
!  (68)  [     0.] Propsi
!  (69)  [     0.] Radstr
!  (70)  [     0.] Nplant
!  (71)  [     0.] Lwavp
! OVHEXP(NMOVAR)   upper expected limit of output quantity
!  ( 1)  [  1.E10] Xp
!  ( 2)  [  1.E10] Yp
!  ( 3)  [  1.E10] Dist
!  ( 4)  [   100.] Depth
!  ( 5)  [     2.] Vel
!  ( 6)  [     1.] Ubot
!  ( 7)  [   100.] Dissip
!  ( 8)  [   100.] Qb
!  ( 9)  [   100.] Leak
!  (10)  [    10.] Hs
!  (11)  [   100.] Tm01
!  (12)  [   100.] Tpeak
!  (13)  [   360.] Dir
!  (14)  [   360.] PkDir
!  (15)  [   360.] TDir
!  (16)  [    60.] Dspr
!  (17)  [   200.] Wlen
!  (18)  [    0.1] Steepn
!  (19)  [    10.] Transp
!  (20)  [    10.] WForce
!  (21)  [   100.] AcDens
!  (22)  [   100.] EnDens
!  (23)  [  1.E10] Aux
!  (24)  [   100.] Xc
!  (25)  [   100.] Yc
!  (26)  [    50.] Windv
!  (27)  [     1.] FrCoef
!  (28)  [   100.] RTm01
!  (29)  [   100.] EnDens
!  (30)  [     1.] dHs
!  (31)  [     2.] dTm
!  (32)  [   100.] Tm02
!  (33)  [     1.] FSpr
!  (34)  [     1.] Urms
!  (35)  [     1.] Ufric
!  (36)  [     1.] Zlen
!  (37)  [     1.] TauW
!  (38)  [     1.] Cdrag
!  (39)  [     1.] Setup
!  (47)  [   100.] Tm_10
!  (48)  [   100.] RTm_10
!  (49)  [    10.] DifPar
!  (50)  [   100.] TmBot
!  (51)  [   100.] Watlev
!  (52)  [   100.] Botlev
!  (53)  [   100.] TPSmoo
!  (54)  [    0.1] Sfric
!  (55)  [    0.1] Ssurf
!  (56)  [    0.1] Swcap
!  (57)  [    0.1] Sveg
!  (58)  [     1.] Qp
!  (59)  [  1000.] BFI
!  (60)  [    0.1] Genera
!  (61)  [    0.1] Swind
!  (62)  [    0.1] Redist
!  (63)  [    0.1] Snl4
!  (64)  [    0.1] Snl3
!  (65)  [    0.1] Propag
!  (66)  [    0.1] Propxy
!  (67)  [    0.1] Propth
!  (68)  [    0.1] Propsi
!  (69)  [    0.1] Radstr
!  (70)  [   100.] Nplant
!  (71)  [   200.] Lwavp
! OVULIM(NMOVAR)   upper limit of validity
!  ( 1)  [  1.E10] Xp
!  ( 2)  [  1.E10] Yp
!  ( 3)  [  1.E10] Dist
!  ( 4)  [   1.E4] Depth
!  ( 5)  [   100.] Vel
!  ( 6)  [    10.] Ubot
!  ( 7)  [  1000.] Dissip
!  ( 8)  [     1.] Qb
!  ( 9)  [  1000.] Leak
!  (10)  [   100.] Hs
!  (11)  [  1000.] Tm01
!  (12)  [  1000.] Tpeak
!  (13)  [   360.] Dir
!  (14)  [   360.] PkDir
!  (15)  [   360.] TDir
!  (16)  [   360.] Dspr
!  (17)  [  1000.] Wlen
!  (18)  [     1.] Steepn
!  (19)  [   100.] Transp
!  (20)  [   1.E5] WForce
!  (21)  [  1000.] AcDens
!  (22)  [  1000.] EnDens
!  (23)  [  1.E10] Aux
!  (24)  [  1000.] Xc
!  (25)  [  1000.] Yc
!  (26)  [   100.] Windv
!  (27)  [     1.] FrCoef
!  (28)  [  1000.] RTm01
!  (29)  [  1000.] EnDens
!  (30)  [   100.] dHs
!  (31)  [   100.] dTm
!  (32)  [  1000.] Tm02
!  (33)  [     1.] FSpr
!  (34)  [    10.] Urms
!  (35)  [    10.] Ufric
!  (36)  [     1.] Zlen
!  (37)  [    10.] TauW
!  (38)  [     1.] Cdrag
!  (39)  [     1.] Setup
!  (47)  [  1000.] Tm_10
!  (48)  [  1000.] RTm_10
!  (49)  [    50.] DifPar
!  (50)  [  1000.] TmBot
!  (51)  [   1.E4] Watlev
!  (52)  [   1.E4] Botlev
!  (53)  [  1000.] TPSmoo
!  (54)  [  1000.] Sfric
!  (55)  [  1000.] Ssurf
!  (56)  [  1000.] Swcap
!  (57)  [  1000.] Sveg
!  (58)  [     1.] Qp
!  (59)  [  1000.] BFI
!  (60)  [  1000.] Genera
!  (61)  [  1000.] Swind
!  (62)  [  1000.] Redist
!  (63)  [  1000.] Snl4
!  (64)  [  1000.] Snl3
!  (65)  [  1000.] Propag
!  (66)  [  1000.] Propxy
!  (67)  [  1000.] Propth
!  (68)  [  1000.] Propsi
!  (69)  [  1000.] Radstr
!  (70)  [  1000.] Nplant
!  (71)  [  1000.] Lwavp
! SINCQ  [       ] sin of ALCQ
! SINPQ  [       ] sin of ALPQ
! SPCPOW [      1] power in expression for computation of average frequency
!                  =power; set by command SET ... [power] ... (not documented)
! XQLEN  [       ] length of x-side of output frame
! XQP    [       ] x-coordinate (user coord.) of origin of output frame
! YQLEN  [       ] length of y-side of output frame
! YQP    [       ] y-coordinate (user coord.) of origin of output frame
!
      INTEGER ERRPTS
      INTEGER INRHOG,         IUBOTR
      INTEGER OVSVTY(NMOVAR), SPCPOW
      REAL    AKPOWR,         ALCQ,        ALPQ,          COSCQ
      REAL    COSPQ,          DXK,         DYK
      REAL    OVEXCV(NMOVAR),              OVLEXP(NMOVAR)
      REAL    OVLLIM(NMOVAR),              OVHEXP(NMOVAR)
      REAL    OVULIM(NMOVAR),              SINCQ,         SINPQ
      REAL    XPQ,            XQLEN,       XQP,           YPQ
      REAL    YQLEN,          YQP,         OUTPAR(MOUTPA)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE SWCOMM1

      MODULE SWCOMM2
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
!     40.41: Marcel Zijlema
!     41.13: Nico Booij
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file SWCOMM2.INC
!     41.13, Jul. 10: NWAMN unused, replaced by LWDATE (length of WAM date)
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     NUMGRD [  11]  maximum number of grids
!
      INTEGER NUMGRD
      PARAMETER (NUMGRD = 11)
!
!  7. Local variables
!
!     *** location and dimensions of input grids ***
!
! The following input grids are known in SWAN:
! grid(1): depth
! grid(2): current velocity x-component
! grid(3): current velocity y-component
! grid(4): friction coefficient
! grid(5): wind velocity x-component
! grid(6): wind velocity y-component
! grid(7): water level
! grid(8): y-coordinate
! grid(9): x-coordinate
! grid(10): air-sea temperature difference
! grid(11): number of plants per square meter
!
! ALPG(NUMGRD)  [    0.] direction of the x-axis w.r.t. user coordinates
! COSPG(NUMGRD) [    1.] cos of ALPG
! COSVC         [CALCUL] =COS(-ALPG(2)),
!                        cos of angle of current input grid w.r.t. user coordinates
! COSWC         [CALCUL] =COS(-ALPG(5)),
!                        cos of angle of wind input grid w.r.t. user coordinates
! CVLEFT        [      ] the (curv.lin.) comput. grid is left/right-oriented

! DXG(NUMGRD)   [    0.] mesh size of input grid in x-direction
! DYG(NUMGRD)   [    0.] mesh size of input grid in y-direction
! DYNDEP        [   FAL] is True if depth varies with time
! EXCFLD(NUMGRD)[-1.E20] exception values for input grids
! ICOND         [     0] initial conditions
!                        =0 when mode stationary, or no initial conditions needed
!                        =1 when mode non-stationary and initial conditions should
!                           be calculated
! IGTYPE(NUMGRD)[      ] =0 when grid has constant values
!                        =1 when grid is regular
!                        =2 when grid is curvilinear
! LEDS(NUMGRD)  [     0] =0 when values have not been read
!                        =1 if values were read
! LWDATE        [     0] length of date-time in WAM file
! LXOFFS        [   FAL] offset values were/were not initialized already
! MXG(NUMGRD)   [     0] number of meshes in x-direction
! MYG(NUMGRD)   [     0] number of meshes in y-direction
! OPTG          [     1] type of the computational grid
!                        =1 when regular
!                        =2 when irregular, but rectangular (not used)
!                        =3 when curvilinear
!                        =5 when unstructured grid
! NBFILS        [     0] number of boundary condition files
! NBSPEC        [     0] number of boundary spectra
! NBGRPT        [     0] number of comp. grid points for which boundary
!                        condition holds
! RDTIM         [      ] =0 when in stationary mode
!                        =1/DT when in non-stationary mode
! SINPG(NUMGRD) [    0.] sin of ALPG
! SINVC         [CALCUL] =SIN(-ALPG(2)),
!                        sin of angle of current input grid w.r.t. comp. grid
! SINWC         [CALCUL] =SIN(-ALPG(5)),
!                        sin of angle of wind input grid w.r.t. comp. grid
! STAGX(NUMGRD) [    0.] staggering of curv.lin. inp. grid w.r.t. comp. grid in X
! STAGY(NUMGRD) [    0.] staggering of curv.lin. inp. grid w.r.t. comp. grid in Y
! VARAST        [   FAL] air-sea temp. diff. is/is not variable over space
! VARFR         [   FAL] friction coefficient is/is not variable over space
! VARNPL        [   FAL] number of plants / m2 is/is not variable over space
! VARWI         [   FAL] wind velocity is/is not variable over space
! VARWLV        [   FAL] water level is/is not variable over space
! XPG(NUMGRD)   [    0.] x of origin
! XOFFS         [    0.] offset value in x
!                        (from user coord. system to internal coord. system)
! YPG(NUMGRD)   [    0.] y of origin
! YOFFS         [    0.] offset value in y
!                        (from user coord. system to internal coord. system)
!
      INTEGER ICOND
      INTEGER IGTYPE(NUMGRD),           LEDS(NUMGRD)
      INTEGER MXG(NUMGRD), MYG(NUMGRD), LWDATE
      INTEGER OPTG,        NBFILS,      NBSPEC,      NBGRPT
      REAL    ALPG(NUMGRD),             COSPG(NUMGRD)
      REAL    COSVC,       COSWC,       DXG(NUMGRD)
      REAL    DYG(NUMGRD), EXCFLD(NUMGRD)
      REAL    RDTIM,       SINPG(NUMGRD)
      REAL    SINVC,       SINWC,       STAGX(NUMGRD)
      REAL    STAGY(NUMGRD),            XOFFS
      REAL    XPG(NUMGRD), YOFFS,       YPG(NUMGRD)
      LOGICAL CVLEFT,      LXOFFS,      VARFR
      LOGICAL VARWI,       VARWLV,      DYNDEP
      LOGICAL VARAST
      LOGICAL VARNPL
!
!     *** variables for input files ***
!
!     IFLBEG     begin time of data on file
!     IFLDYN     if =0: data is stationary, if =1: nonstationary
!     IFLEND     end time of data on file
!     IFLFAC     multiplication factor
!     IFLFRM     format string
!     IFLIDL     lay-out in input file
!     IFLIFM     format identifier
!     IFLINT     time interval of data on file
!     IFLNDF     unit ref number of namelist file
!     IFLNDS     unit ref number of data file
!     IFLNHD     number of heading lines per input field
!     IFLNHF     number of heading lines per file
!     IFLTIM     time of last reading
!
      INTEGER      IFLDYN(NUMGRD)
      INTEGER      IFLIDL(NUMGRD), IFLIFM(NUMGRD), IFLNHF(NUMGRD),
     &             IFLNHD(NUMGRD), IFLNDS(NUMGRD), IFLNDF(NUMGRD)
      REAL         IFLBEG(NUMGRD), IFLINT(NUMGRD), IFLEND(NUMGRD),
     &             IFLTIM(NUMGRD)
      REAL         IFLFAC(NUMGRD)
      CHARACTER*40 IFLFRM(NUMGRD)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE SWCOMM2

      MODULE SWCOMM3
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file SWCOMM3.INC
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     MBOT   [  10] dimension of array PBOT
!     MDIFFR [  10] dimension of array PDIFFR                             40.21
!     MDISP  [   5] dimension for dissipation arrays                      40.67
!     MGENR  [   1] dimension for generation arrays                       40.85
!     MICMAX [  10] max. number of points in comput. stencil
!     MNUMS  [  35] dimension of array PNUMS
!     MQUAD  [  10] dimension of array PQUAD
!     MREDS  [   2] dimension for redistribution arrays                   40.85
!     MSETUP [   2] dimension of array PSETUP
!     MSHAPE [   5] dimension of array PSHAPE
!     MSPPAR [   5] dimension of array SPPARM
!     MSURF  [  15] dimension of array PSURF                              41.06
!     MTRNP  [   3] dimension for propagation arrays                      40.85
!     MTRIAD [  10] dimension of array PTRIAD
!     MWCAP  [  15] dimension of array PWCAP
!     MWIND  [  40] dimension of array PWIND
!
      INTEGER             MBOT,        MNUMS,       MQUAD
      INTEGER             MSETUP,      MSHAPE,      MSPPAR,      MSURF
      INTEGER             MTRIAD,      MWCAP,       MWIND
      INTEGER             MICMAX
      INTEGER             MDIFFR                                          40.21
      INTEGER             MDISP                                           40.67
      INTEGER             MGENR, MREDS, MTRNP                             40.85
!
      PARAMETER           (MBOT   = 10)
      PARAMETER           (MDIFFR = 10)                                   40.21
      PARAMETER           (MDISP  =  5)                                   40.67
      PARAMETER           (MGENR  =  1)                                   40.85
      PARAMETER           (MREDS  =  2)                                   40.85
      PARAMETER           (MTRNP  =  3)                                   40.85
      PARAMETER           (MICMAX = 10)
      PARAMETER           (MNUMS  = 35)
      PARAMETER           (MQUAD  = 10)
      PARAMETER           (MSETUP =  2)
      PARAMETER           (MSHAPE =  5)
      PARAMETER           (MSPPAR =  5)
      PARAMETER           (MSURF  = 15)                                   41.06
      PARAMETER           (MTRIAD = 10)
      PARAMETER           (MWCAP  = 15)
      PARAMETER           (MWIND  = 40)
!
!  7. Local variables
!
!     *** pointers for data arrays on computational grid ***
!
! JABIN  [  1] within array LSWMAT
! JABLK  [  2] within array LSWMAT
! JAOLD  [ 8+2*MDISP] within array SWMATR
! JASTD2 [  1] new air-sea temp. diff. within array COMPDA
! JASTD3 [  1] last read air-sea temp. diff. within array COMPDA
! JBOTLV [  1] bottom level within array COMPDA
!              set by command TABLE/BLOCK
! JCDRAG [  1] drag coefficient within array COMPDA,
!              set by command WCAP JANS ...
! JDHS   [  6] wave height correction within array COMPDA
!              (difference in Hs between last two iterations)
! JDIS0  [  7] within array SWMATR
! JDIS1  [  7+MDISP] within array SWMATR
! JDISS  [  2] total dissipation within array COMPDA
! JDPSAV [  1] saved depth (for setup) within array COMPDA
! JDP1   [  7] old depth within array COMPDA
! JDP2   [  8] new depth within array COMPDA
! JDP3   [ 15] last read depth within array COMPDA
! JDSXB  [  1] bottom friction dissipation within array COMPDA
!              set by command TABLE/BLOCK
! JDSXS  [  1] surf dissipation within array COMPDA
!              set by command TABLE/BLOCK
! JDSXV  [  1] vegetation dissipation within array COMPDA
!              set by command TABLE/BLOCK
! JDSXW  [  1] whitecapping dissipation within array COMPDA
!              set by command TABLE/BLOCK
! JDTM   [ 20] wave period correction within array COMPDA
!              (difference in average wave period between last two iterations)
! JFRC2  [  1] friction coefficient within array COMPDA
!              set by command READ FR ...
! JFRC3  [  1] friction coefficient within array COMPDA
!              set by command READ FR ...
! JGEN0  [   ] within array SWMATR
! JGEN1  [   ] within array SWMATR
! JGENR  [  1] total generation within array COMPDA
!              set by command TABLE/BLOCK
! JGSXW  [  1] wind input within array COMPDA
!              set by command TABLE/BLOCK
! JHS    [  1] significant wave height Hs within array COMPDA
! JHSIBC [ 25] significant wave height from boundary condition in array COMPDA
! JLEAK  [ 21] "leak" within array COMPDA
!              (refractive energy tranport over sector boundaries)
! JLEK1  [  7+2*MDISP] within array SWMATR
! JMAT5  [  5] within array SWMATR
! JMAT6  [  6] within array SWMATR
! JMATD  [  1] within array SWMATR
! JMATL  [  3] within array SWMATR
! JMATR  [  2] within array SWMATR
! JMATU  [  4] within array SWMATR
! JNPLA2 [   ] new number of plants / m2 within array COMPDA
! JNPLA3 [   ] last read number of plants / m2 within array COMPDA
! JP4D   [  7] within array SWTSDA, quadruplet interactions
! JP4S   [  6] within array SWTSDA, quadruplet interactions
! JPBOT  [  1] bottom wave period within array COMPDA
!              set by command TABLE/BLOCK
! JPBTFR [  4] within array SWTSDA, bottom friction
! JPTRI  [  8] within array SWTSDA, triad interactions
! JPVEGT [  9] within array SWTSDA, vegetation dissipation
! JPWBRK [  5] within array SWTSDA, surf breaking
! JPWCAP [  3] within array SWTSDA, white capping
! JPWNDA [  1] within array SWTSDA, wind source term part A
! JPWNDB [  2] within array SWTSDA, wind source term part B
! JQB    [  4] fraction of breaking waves within array COMPDA
! JRADS  [  1] radiation stress within array COMPDA
!              set by command TABLE/BLOCK
! JRED0  [   ] within array SWMATR
! JRED1  [   ] within array SWMATR
! JREDS  [  1] total redistribution within array COMPDA
!              set by command TABLE/BLOCK
! JRSXQ  [  1] quadruplets within array COMPDA
!              set by command TABLE/BLOCK
! JRSXT  [  1] triads within array COMPDA
!              set by command TABLE/BLOCK
! JSETUP [  1] setup values within array COMPDA
! JSTP   [  5] steepness within array COMPDA
! JTAUW  [  1] TauW within array COMPDA,
!              set by command WCAP JANS ...
! JTRA0  [   ] within array SWMATR
! JTRA1  [   ] within array SWMATR
! JTRAN  [  1] total propagation within array COMPDA
!              set by command TABLE/BLOCK
! JTSXG  [  1] xy-propagation within array COMPDA
!              set by command TABLE/BLOCK
! JTSXT  [  1] theta-propagation within array COMPDA
!              set by command TABLE/BLOCK
! JTSXS  [  1] sigma-propagation within array COMPDA
!              set by command TABLE/BLOCK
! JUBOT  [  3] bottom orbital velocity within array COMPDA
! JURSEL [ 27] Ursell number as used in Triad computation
! JUSTAR [  1] friction velocity within array COMPDA,
!              set by command WCAP JANS ...
! JVX1   [  9] x of old current velocity within array COMPDA
! JVX2   [ 11] x of new current velocity within array COMPDA
! JVX3   [ 13] x of last read current velocity within array COMPDA
! JVY1   [ 10] y of old current velocity within array COMPDA
! JVY2   [ 12] y of new current velocity within array COMPDA
! JVY3   [ 14] y of last read current velocity within array COMPDA
! JWLV1  [ 22] old water level within array COMPDA
! JWLV2  [ 24] new water level within array COMPDA
! JWLV3  [ 23] last read water level within array COMPDA
! JWX2   [ 16] x of new wind velocity within array COMPDA
! JWX3   [ 18] x of last read wind velocity within array COMPDA
! JWY2   [ 17] y of new wind velocity within array COMPDA
! JWY3   [ 19] y of last read wind velocity within array COMPDA
! JZEL   [  1] roughness within array COMPDA,
!              set by command WCAP JANS ...
! MCMVAR [ 27] within array COMPDA,
!              =MCMVAR+2, for command READ FR ... (add JFRC2, JFRC3)
!              =MCMVAR+4, for command WCAP JANS ... (add JCDRAG, JTAUW,
!               JUSTAR, JZEL)
!              =MCMVAR+1, for command TABLE/BLOCK TMBOT (add JPBOT)
!              =MCMVAR+1, for command TABLE/BLOCK BOTLEV (add JBOTLV)
!              =MCMVAR+1, for command TABLE/BLOCK DISB (add JDSXB)
!              =MCMVAR+1, for command TABLE/BLOCK DISSU (add JDSXS)
!              =MCMVAR+1, for command TABLE/BLOCK DISW (add JDSXW)
!              =MCMVAR+1, for command TABLE/BLOCK DISV (add JDSXV)
!              =MCMVAR+1, for command TABLE/BLOCK GENE (add JGENR)
!              =MCMVAR+1, for command TABLE/BLOCK GENW (add JGSXW)
!              =MCMVAR+1, for command TABLE/BLOCK REDI (add JREDS)
!              =MCMVAR+1, for command TABLE/BLOCK REDQ (add JRSXQ)
!              =MCMVAR+1, for command TABLE/BLOCK REDT (add JRSXT)
!              =MCMVAR+1, for command TABLE/BLOCK PROPA (add JTRAN)
!              =MCMVAR+1, for command TABLE/BLOCK PROPX (add JTSXG)
!              =MCMVAR+1, for command TABLE/BLOCK PROPT (add JTSXT)
!              =MCMVAR+1, for command TABLE/BLOCK PROPS (add JTSXS)
!              =MCMVAR+1, for command TABLE/BLOCK RADST (add JRADS)
! MLSWMAT [  2] within array LSWMAT
! MSWMATR [ 14] within array SWMATR
! MTSVAR [  9] within array TESTDA
!
      INTEGER             JABIN,       JABLK,       JAOLD
      INTEGER             JASTD2,      JASTD3
      INTEGER             JCDRAG,      JDHS,        JDIS0
      INTEGER             JDIS1,       JDISS,       JDPSAV,      JDP1
      INTEGER             JDP2
      INTEGER             JDP3,        JDTM,        JFRC2,       JFRC3
      INTEGER             JGEN0,       JGEN1,       JRED0,       JRED1
      INTEGER             JTRA0,       JTRA1
      INTEGER             JDSXB
      INTEGER             JDSXS
      INTEGER             JDSXW
      INTEGER             JDSXV
      INTEGER             JGSXW,       JGENR
      INTEGER             JRSXQ,       JRSXT,       JREDS,       JRADS
      INTEGER             JTSXG,       JTSXT,       JTSXS,       JTRAN
      INTEGER             JHS,         JHSIBC,      JLEAK
      INTEGER             JLEK1,       JMAT5,       JMAT6,       JMATD
      INTEGER             JMATL,       JMATR,       JMATU
      INTEGER             JP4D,        JP4S,        JPBTFR,      JPTRI
      INTEGER             JPWBRK,      JPWCAP,      JPWNDA,      JPWNDB
      INTEGER             JPVEGT,      JNPLA2,      JNPLA3
      INTEGER             JQB,         JSETUP,      JSTP,        JTAUW
      INTEGER             JUBOT,       JUSTAR,      JVX1,        JVX2
      INTEGER             JVX3,        JVY1,        JVY2,        JVY3
      INTEGER             JWLV1,       JWLV2,       JWLV3
      INTEGER             JWX2,        JWX3
      INTEGER             JWY2,        JWY3,        JZEL  ,      JPBOT
      INTEGER             MCMVAR,                   MTSVAR,      JURSEL
      INTEGER             JBOTLV,      MSWMATR,     MLSWMAT
!
!     *** location and dimensions of computational grid ***
!
! ALCP   [CALCUL] =-ALPC;
!                  direction of user coordinates w.r.t. computational coordinates
! ALOCMP [ FALSE]  if True, array COMPDA must be re-allocated
! ALPC   [    0.]  direction of x-axis of computational grid w.r.t. user coordinates
! COSPC  [CALCUL] =COS(ALPC)
! DX     [CALCUL] =XCLEN/MXS; mesh size in x-direction of computational grid
! DY     [CALCUL] =YCLEN/MYS; mesh size in y-direction of computational grid
! DDIR   [CALCUL] =(SPDIR2-SPDIR1)/MDC;
!                  mesh size in theta-direction of computational grid
! DXRP   [      ]  unused
! DYRP   [      ]  unused
! FRINTF [CALCUL] =ALOG(SHIG/SLOW)/(MSC-1); frequency integration factor (df/f)
!                  (integral over frequency of G(f) = SUM_j (sigma_j*Gj*FRINTF) )
! FRINTH [CALCUL] =SQRT(SFAC); frequency mesh boundary factor
!                  (mesh in frequency space runs from sigma/FRINTH to sigma*FRINTH)
! FULCIR [  TRUE]  spectral directions cover the full/part of circle
! ICOMP  [     1]  unused
! ILMAX  [CALCUL]  maximum number of layers to be used in vegetation model
! IXCGRD [      ]  IX of points of computational stencil
! IYCGRD [      ]  IY of points of computational stencil
! KCGRD  [      ]  grid address of points of computational stencil
! MCGRD  [     1]  number of wet grid points of the computational grid
! MDC    [     0]  grid points in theta-direction of computational grid
! MDC4MA [CALCUL] =IDHGH; some counter for quadruplet interactions. Stored in WWINT(18)
! MDC4MI [CALCUL] =IDLOW; some counter for quadruplet interactions. Stored in WWINT(17)
! MMCGR  [CALCUL] =MXC*MYC; grid points in computational grid
! MSC    [     0]  grid points in sigma-direction of computational grid
! MSC4MA [CALCUL] =ISHGH; some counter for quadruplet interactions. Stored in WWINT(16)
! MSC4MI [CALCUL] =ISLOW; some counter for quadruplet interactions. Stored in WWINT(15)
! MTC    [     1]  computational time steps
! MXC    [     0]  grid points in x-direction of computational grid
! MYC    [     0]  grid points in y-direction of computational grid
! NGRBND [CALCUL]  number of grid points on computational grid boundary
! NX     [CALCUL] =MXC-1; only used locally. Equal to MXS
! NY     [CALCUL] =MYC-1; only used locally. Equal to MYS
! SHIG   [CALCUL] =2*PI*FRHIG; highest spectral value of sigma
! SINPC  [CALCUL] =SIN(ALPC)
! SLOW   [CALCUL] =2*PI*FRLOW; lowest spectral value of sigma
! SPDIR1 [    0.]  represents first spectral direction (if FULCIR=.FALSE.)
! SPDIR2 [      ]  represents second spectral direction (if FULCIR=.FALSE.)
! XCGMAX [CALCUL]  maximum x-coordinate of computational grid points
! XCGMIN [CALCUL]  minimum x-coordinate of computational grid points
! XCLEN  [      ]  length of computational grid in x-direction
! XCP    [CALCUL] =-XPC*COSPC-YPC*SINPC;
!                  origin of user coordinates w.r.t. computational coordinates
! XPC    [    0.]  x coordinate of origin of computational grid
! YCGMAX [CALCUL]  maximum y-coordinate of computational grid points
! YCGMIN [CALCUL]  minimum y-coordinate of computational grid points
! YCLEN  [      ]  length of computational grid in y-direction
! YCP    [CALCUL] =XPC*SINPC-YPC*COSPC;
!                  origin of user coordinates w.r.t. computational coordinates
! YPC    [    0.]  y coordinate of origin of computational grid
!
      INTEGER             IXCGRD(MICMAX), IYCGRD(MICMAX), KCGRD(MICMAX)
      INTEGER             ICOMP,       MCGRD
      INTEGER             MDC,         MDC4MA,      MDC4MI,      MMCGR
      INTEGER             MSC,         MSC4MA,      MSC4MI,      MTC
      INTEGER             MXC,         MYC,         NX,          NY
      INTEGER             NGRBND
      INTEGER             ILMAX
      REAL                ALCP,        ALPC,        COSPC,       DDIR
      REAL                DX,          DXRP,        DY,          DYRP
      REAL                FRINTF,      FRINTH,      SHIG,        SINPC
      REAL                SLOW,        SPDIR1,      SPDIR2,      XCLEN
      REAL                XCP,         XPC,         YCLEN,       YCP
      REAL                YPC
      REAL                XCGMIN,      XCGMAX,      YCGMIN,      YCGMAX
      LOGICAL             FULCIR
      LOGICAL ::          ALOCMP = .FALSE.                                40.97
!$OMP THREADPRIVATE(IXCGRD,IYCGRD,KCGRD)
!
!     *** physical parameters ***
!
! CASTD  [    0.] (constant) air-sea temperature difference
! CDCAP  [99999.] maximum drag coefficient
! DEGRAD [CALCUL] =PI/180; constant to transform degrees to radians
! DNORTH [   90.] direction of the North w.r.t. x-axis of user coordinates
!                 =nor; set by command SET ... [nor] ...
! GRAV   [  9.81] acceleration due to gravity
!                 =grav; set by command SET ... [grav] ...
! PI     [3.1415] circular constant
! PI2    [CALCUL] =2*PI;
! PWTAIL(10)      coefficients to calculate tail of the spectrum
!    (1) [    4.] tail power of energy density spectrum as function of freq.
!                 =5.; for command GEN3 JANS ...
!                 =5.; for command WCAP JANS ...,
!                      not documented in manual
!                 =5.; for command GROWTH G3 JANS ...,
!                      not documented in manual
!                 =pwtail; set by command SET ... [pwtail]
!    (2) [   2.5] energy spectrum w.r.t. wave number, not used
!    (3) [CALCUL] =PWTAIL(1)+1; action density spectrum w.r.t. frequency
!    (4) [    3.] action density spectrum w.r.t. wave number, not used
!    (5) [CALCUL] =1./(PWTAIL(1)*(1.+PWTAIL(1)*(FRINTH-1.)));
!                 tail factor for the action integral
!    (6) [CALCUL] =1./((PWTAIL(1)-1.)*(1.+(PWTAIL(1)-1.)*(FRINTH-1.)));
!                 tail factor for the energy integral
!    (7) [CALCUL] =1./((PWTAIL(1)-2.)*(1.+(PWTAIL(1)-2.)*(FRINTH-1.)));
!                 tail factor for the first moment of energy
!    (8) [CALCUL] =1./((PWTAIL(1)-3.)*(1.+(PWTAIL(1)-3.)*(FRINTH-1.)));
!                 tail factor for the second moment of energy
! RHO    [ 1025.] density of the water
!                 =rho; set by command SET ... [rho] ...
! WLEV   [    0.] water level
!                 =level; set by command SET [level] ...
!
      REAL CASTD,       CDCAP
      REAL DEGRAD,      DNORTH,      GRAV,        PI
      REAL PI2,         PWTAIL(10),  RHO,         WLEV
!
!     *** information related to the numerical scheme ***
!
! ACUPDA [  true] indicates whether or not action densities are to be updated
!                 during computation
! BNAUT  [ false] indicates whether nautical or cartesian directions are used
! BNDCHK [  true] indicates whether computed Hs on boundary must be compared
!                 with value entered as boundary condition
! BRESCL [  true] rescaling on/off
! CSETUP [  true] indicates whether solver for setup has converged
! DEPMIN [  0.05] threshold depth (to prevent zero divisions)
!                 =depmin; set by command SET ... [depmin] ...
! DSHAPE [     2] indicates option for computation of directional distribution
!                 in the spectrum (boundary spectra etc.)
!                 =1: directional spread in degrees is given
!                 =2: power of COS is given
! FSHAPE [     2] indicates option for computation of frequency distribution
!                 in the spectrum (boundary spectra etc.)
!                 =2: Jonswap(default set by subr SWINIT),
!                 =1: Pierson-Moskowitz, =3: bin, =4: Gaussian (set by command
!                 BOUNshape ..)
! HSRERR [   0.1] The error margin allowed between pre-scribed and calculated Hs
!                 at the up-wave boundary. If exceeded a warning is produced
! IBOT   [     0] indicator bottom friction:
!                 =0; no bottom friction dissipation
!                 =1; set by command FRIC JON  ..., Jonswap bottom friction model
!                 =2; for command FRIC COLL ..., Collins bottom friction model
!                 =3; for command FRIC MAD  ..., Madsen bottom friction model
! ICMAX  [     3] number of points in computational stencil
! ICOR   [      ] not used
! ICUR   [     0] indicates presence of currents:
!                 =0; no currents
!                 =1; for command READ CUR ..., currents are present
! IDBR   [     1] not used
! IDIF   [     0] not used
! IDIFFR [     0] diffraction method                                      40.21
!                 0= no diffraction                                       40.21
!                 1= diffraction                                          40.21
! IFRSRF [     0] indicates frequency dependent surf breaking             41.06
!                 0=no frequency dependency                               41.06
!                 1=frequency dependency                                  41.06
! IGEN   [     3] indicates the generation mode
!                 =1; for command GEN1 ...,
!                 =2; for command GEN2 ...,
!                 =3; for command GEN3 ...,
! IINC   [     0] not used
! IPRE   [      ] not used
! IQUAD  [     2] indicates quadruplet interaction term:
!                 =0; for command OFF QUAD
!                 =0; for command GEN1 ...,
!                 =0; for command GEN2 ...,
!                 =0; for command GROWTH G1 ... (not documented in manual),
!                 =0; for command GROWTH G2 ... (not documented in manual),
!                 quadruplets are inactive
!                 =1; quadruplets are calculated semi implicit per sweep direction
!                 =2; for command GEN3 ...,
!                 =2; for command QUAD, not documented in the manual
!                 =2; set when IWIND=3 or 4 and ICUR=0 in SUBR ERRCHK,
!                 quadruplets are calculated fully explicit per sweep direction
!                 =3; set when IWIND=3 or 4 and ICUR=1 in SUBR ERRCHK,
!                 quadruplets are calculated fully explicit per iteration
!                 =8; quadruplets are calculated fully explicit per iteration and
!                     interactions are interpolated in piecewise constant manner
!                 =iquad; set by command GEN3 ... QUAD [iquad] ...,
! IREFR  [     1] indicates refraction effect:
!                 =0; for command OFF REF, refraction is inactive
!                 =1; refraction is active
! ISURF  [     1] indicates surf breaking (shallow water) term:
!                 =0; for command OFF BRE, surf breaking is inactive
!                 =1; for command BRE CON ..., surf breaking with constant parameter
!                 =2; for command BRE VAR ..., surf breaking according to Nelson
!                 =3; for command BRE RUE ..., surf breaking according to Ruessink
!                 =4; for command BRE TG ... , surf breaking according to Thornton and Guza
! ITERMX [      ] maximum number of iterations:
!                 is set equal to MXITST in case of stationary computations
!                 is set equal to MXITNS in case of nonstationary computations
! ITFRE  [     1] indicator for transport of action in frequency space
!                 =0; for command OFF FSH, frequency shifting inactive
!                 =1; frequency shifting active
! ITRIAD [     0] indicates triad interaction term:
!                 =0; triads are inactive
!                 =1; for command TRI DTA IMP ..., not documented in manual
!                 =2; for command TRI DTA EXP ..., not documented in manual
!                 =3; for command TRI [trfac] [cutfr], as in manual
!                 =3; for command TRI LTA IMP ..., not documented in manual
!                 =4; for command TRI LTA EXP ..., not documented in manual
! IVEG   [     0] indicates vegetation dissipation term:
!                 =0; no dissipation due to vegetation
!                 =1; dissipation due to vegetation according to Dalrymple (1984)
! IWCAP  [     1] indicates whitecapping:
!                 =0; for command GEN1 ...,
!                 =0; for command GEN2 ...,
!                 =0; for command OFF WCAP, no whitecapping
!                 =1; for command GEN3 KOM ...,
!                 =1; for command WCAP KOM ..., not documented in manual,
!                 standard WAM formulation (Komen et al.; 1984)
!                 =2; for command GEN3 JANS ...,
!                 =2; for command WCAP JANS ..., not documented in manual,
!                 according to Janssen (1989, 1991)
!                 =3; for command WCAP LHIG ..., not documented in manual,
!                 according to Longuet-Higgins (1967), Yuan et al. (1986)
!                 =4; for command WCAP BJ ..., not documented in manual,
!                 according to Battjes & Janssen (1978)
!                 =5; for command WCAP KBJ ..., not documented in manual,
!                 combined formulation of Komen (1) and Battjes & Janssen (4)
! IWIND  [     0] indicates presence of wind, and type of source term used:
!                 =0; no wind
!                 =1; for command GEN1 ..., if wind is made active,
!                 =1; for command GROWTH G1 ..., not documented in manual,
!                 1st generation source term
!                 =2; for command GEN2 ..., if wind is made active,
!                 =2; for command GROWTH G2 ..., not documented in manual,
!                 2nd generation source term (as in Dolphin)
!                 =3; for command WIND ..., if IWIND still was 0, else unchanged,
!                 =3; for command GEN3 KOM ..., if wind is made active,
!                 =3; for command GROWTH G3 KOM ..., not documented in manual,
!                 3rd generation source term (Snyder)
!                 =4; for command GEN3 JANS ..., if wind is made active,
!                 =4; for command GROWTH G3 JANS ..., not documented in manual,
!                 source term by P. Janssen (1989, 1991)
!                 =5; for command GEN3 YAN ..., if wind is made active,
!                 =5; for command GROWTH G3 YAN ..., not documented in manual,
!                 not documented in manual
! LADDS  [.fals.] indicates whether extra output is requested or not
! LSETUP [     0] =0; setup is not calculated
!                 =1; setup is calculated
!                 =2; setup is calculated with the boundary conditions from
!                     a nest file
! MXITST [    50] max. number of iterations in stationary computations
! MXITNS [     1] max. number of iterations in nonstationary computations
! NCOR   [     1] not used
! NSTATC [     1] indicates stationarity of computation:
!                 =0; stationary computation
!                 =1; nonstationary computation
! NSTATM  [   -1] 0: stationary mode, 1: nonstationary mode, -1: unknown
! NUMOBS [      ] number of obstacles
! NCOMPT [      ] number of COMPUTE commands
! OFFSRC [.fals.] indicates whether source terms are included in
!                 action balance equation or not
! ONED   [.FALS.] Indicates whether the calculation should be performed in 1D-mode
! PDIFFR(MDIFFR)  coefficients for diffraction                            40.21
!  ( 1)           smoothing parameter                                     40.21
!  ( 2)           number of smoothing steps                               40.21
!  ( 3)           if (Cx,Cy)-velocities are modified or not               40.21
! PBOT(MBOT)      coefficients for the bottom friction models
!  ( 1)  [    0.] =cfc; set by command FRIC COL [cfw] [cfc], (Collins equation),
!                 not documented in the manual
!  ( 2)  [ 0.015] =cfw; set by command FRIC COL [cfw], (Collins equation),
!                 also used as CFW in the source code
!  ( 3)  [ 0.067] =cfjon; set by command FRIC JON [cfjon], (Jonswap equation)
!  ( 4)  [ -0.08] =mf; value cannot be changed, (Madsen equation)
!  ( 5)  [  0.05] bottom roughness length scale, (Madsen equation),
!                 =kn; set by command FRIC JON [kn],
!                 also used as AKN in the source code
! PNUMS(MNUMS)    numerical coefficients
!                 accuracy criterion:
!  ( 1)  [  0.02] relative error in Hs and Tm01
!                 =drel; set by command NUM ACCUR [drel] ...
!  ( 2)  [  0.03] absolute error in Hs (m)
!                 =dhabs; set by command NUM ACCUR ... [dhabs] ...
!  ( 3)  [   0.3] absolute error in Tm01 (s)
!                 =dtabs; set by command NUM ACCUR ... [dtabs] ...
!  ( 4)  [ 98.00] percentage of wet grid points were absolute and relative
!                 accuracy has been reached
!                 =npnts; set by command NUM ACCUR ... [npnts] ...
!  ( 5)  [    0.] not used
!                 diffusion schemes:
!  ( 6)  [   0.5] numerical diffusion over theta,
!                 =cdd; set by command NUM DIR [cdd]
!  ( 7)  [   0.5] numerical diffusion over sigma,
!                 =css; set by command NUM SIGIM [css] ...
!  ( 8)  [    1.] numerical scheme in frequency space
!                 =1.; for command NUM SIGIM, implicit scheme
!                 =2.; for command NUM SIGEX, explicit scheme,
!                 with CFL criterion
!                 =3.; for command NUM FIL, explicit scheme,
!                 without CFL criterion, not documented in the manual
!  ( 9)  [  0.01] diffusion coefficient for explicit scheme,
!                 =diffc; set by command NUM FIL [diffc],
!                 not documented in the manual
!                 iterative solver (SIP):
!  (12)  [ 1.E-4] termination criterion for iterative solver
!                 =eps2; set by command NUM SIGIM ... [eps2] ...
!                 (||Ax-b|| < eps2 * ||b||)
!  (13)  [    0.] output for SIP solver
!                 =outp; set by command NUM SIGIM ... [outp] ...
!                 <0. no output
!                 =0. only fatal errors are printed
!                 =1. additional information about the iteration is printed
!                 =2. maximal output regarding the iteration process
!  (14)  [   20.] maximum number of iterations in the solver
!                 =niter; set by command NUM SIGIM ... [niter]
!  (15)  [  0.02] global error in Hs
!  (16)  [  0.02] global error in Tm01
!  (17)  [   -1.] coefficient for limitation of Ctheta (not used currently)
!  (18)  [   0.8] limitation on Froude number (current velocity is reduced if
!                 larger than CGMAX=PNUMS(18)*SQRT(GRAV*DEPW)),
!                 =froudmax; set by command SET ... [froudmax] ...,
!                 not documented in the manual
!  (19)  [CALCUL] =0.5*SQRT(2.); CFL criterion for explicit scheme in
!                 frequency space, manual mentions default = 0.7,
!                 =cfl; set by command NUM SIGEX [cfl]
!  (20)  [   0.1] maximum growth in spectral bin,
!                 =0.1; for command WIND ..., if IWIND=3 or 4
!                 =0.1; for command GEN3
!                 =0.1; for command TRI
!                 =0.1; for command QUAD,
!                 not documented in the manual
!                 =1.E20; for command OFF QUAD
!                 =limiter; set by command GEN3 ... QUAD ... [limiter]
!                 =limiter; set by command NUM ACCUR ... [itermax] [limiter],
!                 not documented in the manual
!                 =limiter; set by command QUAD [iquad] [limiter],
!                 not documented in the manual
!  (21)  [    0.] =coefficient for type stopping criterion
!  (23)           termination criterion for iterative solver in set-up calculation
!                 =eps2, set by command NUM SETUP ... [eps2] ...
!  (24)           output for iterative solver in set-up calculation
!                 =outp, set by command NUM SETUP ... [outp] ...
!                 <0. no output
!                 =0. only fatal errors are printed
!                 =1. additional information about the iteration is printed
!                 =2. maximal output regarding the iteration process
!  (25)           maximum number of iterations for solver in set-up calculation
!                 =niter, set by command NUM SETUP ... [niter] ...
!  (28)  [    1.] Qb-value at which the limiter is not active in case lowering AC2
!  (30)  [    0.] =under-relaxation factor
! PQUAD(MQUAD)    coefficients for quadruplet interaction
!  ( 1)  [  0.25] lambda in eq. B29 of user manual
!  ( 2)  [  3.E7] coefficient of interactions
!  ( 3)  [   5.5] coefficient for shallow water interactions
!  ( 4)  [ 0.833] coefficient for shallow water interactions
!  ( 5)  [ -1.25] coefficient for shallow water interactions
! PSETUP(MSETUP)
!  ( 1)  [   0.0] not used
!  ( 2)  [   0.0] user defined level for correction of the setup
! PSHAPE(MSHAPE)  coefficients for calculation of spectrum from integral
!                 parameters
!  ( 1)  [   3.3] peak enhancement factor of Jonswap spectrum
!                  =gamma; set by command BOUN SHAPE .. JON [gamma]
!  ( 2)  [   0.1] width of Gaussian spectrum (in Hz)
!                 =sigfr; set by command BOUN SHAPE .. GAU [sigfr]
!                 after reading the value is converted to radians/second
!                 by a 2 pi multiplication
! PSURF(MSURF)    surf breaking coefficients
!  ( 1)  [   1.0] coef. for determining rate of dissipation, (Battjes Janssen),
!                 =1.5; for command BRE VAR
!                 =alpha, set by command BRE CON [alpha] ...
!                 =alpha, set by command BRE VAR [alpha] ...
!  ( 2)  [  0.73] breaker parameter
!                 =gamma, set by command BRE CON ... [gamma]
!  ( 4)  [      ] the min. value of the breaker parameter of Nelson
!                 =0.55; set by command BRE VAR
!                 =gammin; set by command BRE VAR ... [gammin] ...
!  ( 5)  [      ] the max. value of the breaker parameter of Nelson
!                 =0.81; set by command BRE VAR
!                 =gammax; set by command BRE VAR ... [gammax] ...
!  ( 6)  [      ] breaker parameter for negative bottom slopes
!                 =0.73; set by command BRE VAR
!                 =gamneg; set by command BRE VAR ... [gamneg] ...
!  ( 7)  [      ] proportionality coefficient in expression of Nelson
!                 =0.88; set by command BRE VAR
!                 =coeff1; set by command BRE VAR ... [coeff1] ...
!  ( 8)  [      ] coefficient in the exp in the expression of Nelson
!                 =0.012; set by command BRE VAR
!                 =coeff2; set by command BRE VAR ... [coeff2]
! PTRIAD(MTRIAD)
!  ( 1)  [  0.05] controls the proportionality coefficient,
!                 =trfac; set by command TRIAD ... [trfac] ...
!  ( 2)  [   2.5] controls the maximum frequency considered in the comp.,
!                 =cutfr; set by command TRIAD ... [cutfr]
!  ( 3)  [  10.0] controls above which Ursell number the quadruplets (and
!                 thus the limiter) are switched off
!                 =ursell; set by command LIM ... [ursell]
!  ( 4)  [   0.2] critical Ursell number appearing in the biphase expression
!                 =urcrit; set by command TRIAD ... [urcrit]
!  ( 5)  [  0.01] controls below which Ursell number the triads are switched off
!                 =urslim; set by command TRIAD ... [urslim]
! PWCAP(MWCAP)     whitecapping coefficients
!  ( 1)  [2.36E-5] coefficient for Komen et al. (1984),
!                 ALFAWC (Emperical coefficient)
!                 =cds2; set by command GEN3 KOM [cds2] ...,
!                =cds2; set by command WCAP KOM [cds2] ...,
!                 not documented in the user manual
!  ( 2)  [3.02E-3] coefficient for Komen et al. (1984),
!                 ALFAPM (Alpha of Pierson Moskowitz frequency)
!                 =stpm; set by command GEN3 KOM ... [stpm],
!                 =stpm; set by command WCAP KOM ... [stpm],
!                 not documented in the user manual
!  ( 3)  [   4.5] coeff. for Janssen (1989,1991), acc. to Komen et al. (1994),
!                 CFJANS (cds coefficient)
!                 =cds1; set by command GEN3 JANS [cds1] ...,
!                 =cds1; set by command WCAP JANS [cds1] ...,
!                 not documented in the user manual
!  ( 4)  [   0.5] (=DELTA)
!                 =delta; set by command GEN3 JANS ... [delta],
!                 =delta; set by command WCAP JANS ... [delta],
!                 not documented in the user manual
!  ( 5)  [    1.] coefficient of Longuet Higgins,
!                 =cflhig; set by command WCAP LHIG [cflhig],
!                 not documented in the user manual
!  ( 6)  [  0.88] GAMBTJ (Steepness limited wave breaking)
!                 =bjstp; set by command WCAP BJ [bjstp] ...
!                 =bjstp; set by command WCAP KBJ [bjstp] ...
!  ( 7)  [    1.] Alpha in Battjes/Janssen,
!                 =bjalf; set by command WCAP BJ ... [bjalf]
!                 =bjalf; set by command WCAP KBJ ... [bjalf] ...
!  ( 8)  [  0.75] numerical diffusion over sigma
!                 =kconv; set by command WCAP KBJ ... [kconv]
! PWIND(MWIND)    wind growth term coefficients
!  ( 1)  [  188.] controls linear wave growth,
!                 =cf10; set by command GEN1 [cf10] ...
!                 =cf10; set by command GEN2 [cf10] ...
!                 =cf10; set by command GROWTH G1 [cf10] ...,
!                 not documented in the user manual
!                 =cf10; set by command GROWTH G2 [cf10] ...,
!                 not documented in the user manual
!  ( 2)  [  0.59] controls the exponential wave growth,
!                 =cf20; set by command GEN1 ... [cf20] ...
!                 =cf20; set by command GEN2 ... [cf20] ...
!                 =cf20; set by command GROWTH G1 ... [cf20] ...,
!                 not documented in the user manual
!                 =cf20; set by command GROWTH G2 ... [cf20] ...,
!                 not documented in the user manual
!  ( 3)  [  0.12] controls the exponential wave growth,
!                 =cf30; set by command GEN1 ... [cf30] ...
!                 =cf30; set by command GEN2 ... [cf30] ...
!                 =cf30; set by command GROWTH G1 ... [cf30] ...,
!                 not documented in the user manual
!                 =cf30; set by command GROWTH G2 ... [cf30] ...,
!                 not documented in the user manual
!  ( 4)  [  250.] controls the dissipation rate,
!                 =cf40; set by command GEN1 ... [cf40] ...
!                 =cf40; set by command GEN2 ... [cf40] ...
!                 =cf40; set by command GROWTH G1 ... [cf40] ...,
!                 not documented in the user manual
!                 =cf40; set by command GROWTH G2 ... [cf40] ...,
!                 not documented in the user manual
!  ( 5)  [0.0023] controls the spectral energy of the limit spectrum
!                 =cf50; set by command GEN2 ... [cf50] ...
!                 =cf50; set by command GROWTH G2 ... [cf50] ...,
!                 not documented in the user manual
!  ( 6)  [-0.223] controls the spectral energy of the limit spectrum
!                 =cf60; set by command GEN2 ... [cf60] ...
!                 =cf60; set by command GROWTH G2 ... [cf60] ...,
!                 not documented in the user manual
!  ( 7)  [    0.] cf70, not used
!  ( 8)  [ -0.56] cf80, not used
!  ( 9)  [CALCUL] density air / density water (=RHOAW),
!                 =PWIND(16)/RHO
!                 =rhoaw; set by command GROWTH G1 ... [rhoaw] ...,
!                 not documented in the user manual
!                 =rhoaw; set by command GROWTH G2 ... [rhoaw] ...,
!                 not documented in the user manual
!  (10)  [0.0036] limit energy Pierson Moskowitz spectrum
!                 =edmlpm; set by command GEN1 ... [edmlpm] ...
!                 =edmlpm; set by command GEN2 ... [edmlpm] ...
!                 =edmlpm; set by command GROWTH G1 ... [edmlpm] ...,
!                 not documented in the user manual
!                 =edmlpm; set by command GROWTH G2 ... [edmlpm] ...,
!                 not documented in the user manual
!  (11)  [0.00123] drag coefficient
!                 =cdrag; set by command GEN1 ... [cdrag] ...
!                 =cdrag; set by command GEN2 ... [cdrag] ...
!                 =cdrag; set by command GROWTH G1 ... [cdrag] ...,
!                 not documented in the user manual
!                 =cdrag; set by command GROWTH G2 ... [cdrag] ...,
!                 not documented in the user manual
!  (12)  [   1.0] minimum wind velocity, relative to current at 10 m above msl
!                 =umin; set by command GEN1 ... [umin] ...
!                 =umin; set by command GEN2 ... [umin] ...
!                 =umin; set by command GROWTH G1 ... [umin] ...,
!                 not documented in the user manual
!                 =umin; set by command GROWTH G2 ... [umin] ...,
!                 not documented in the user manual
!  (13)  [  0.13] coefficient that determines the Pierson Moskowitz spectrum
!                 =cfpm; set by command GEN1 ... [cfpm]
!                 =cfpm; set by command GEN2 ... [cfpm]
!                 =cfpm; set by command GROWTH G1 ... [cfpm],
!                 not documented in the user manual
!                 =cfpm; set by command GROWTH G2 ... [cfpm],
!                 not documented in the user manual
!  (14)  [  0.01] (=ALPHA) Alpha, according to Janssen (1991) wave growth model
!  (15)  [  0.41] (=XKAPPA)carnock: Kappa
!  (16)  [  1.28] density of the air (=RHOA)
!  (17)  [CALCUL] =RHO, density of the water (=RHOW)
!  (31)  [    0.] proportionality coefficient in the wave growth term of
!                 Caveleri and Malanotte,
!                 =0.0015; for command GEN3 ... AGROW
!                 =a; set by command GEN3 ... AGROW [a]
! SIGMAG [   0.1] width of the Gaussian frequency spectrum in Hz
!                 =0.01; for command BOU STAT ... GAU
!                 =sigfr; set by command BOU STAT ... GAU [sigfr]
!                 after reading the value is converted to radians/second
!                 by a 2 pi multiplication
! SPPARM          integral parameters used for computation of incident spectrum
!  ( 1)  [   ---] significant wave height
!  ( 2)  [   ---] wave period (peak or mean)
!  ( 3)  [   ---] average wave direction
!  ( 4)  [   ---] directional distribution coefficient
! SY0    [   3.3] peak enhancement parameter of the JONSWAP spectrum,
!                 =gamma; set by command BOU STAT ... JON [gamma]
! U10    [    0.] wind velocity
!                 =vel; set by command WIND [vel] ...
! WDIC   [CALCUL] =PI2*((WDIP/PI2-NINT(WDIP/PI2)),
! WDIP   [    0.] wind direction with respect to problem coordinates
!                 =dir; set by command WIND ... [dir]
!
      INTEGER             DSHAPE,      FSHAPE
      INTEGER             ICMAX
      INTEGER             IBOT,        ICOR,        ICUR
      INTEGER             IDBR,        IDIF,        IGEN,        IINC
      INTEGER             IPRE,        IQUAD,       IREFR,       ISURF
      INTEGER             ITERMX,      ITFRE,       ITRIAD
      INTEGER             IWCAP,       IWIND,       LSETUP
      INTEGER             MXITST,      MXITNS,                   NCOR
      INTEGER             NSTATC,      NSTATM,      NUMOBS,      NCOMPT   40.41
      INTEGER             IDIFFR                                          40.21
      INTEGER             IVEG                                            40.55
      INTEGER             IFRSRF                                          41.06
      REAL                DEPMIN,      PBOT(MBOT),  PNUMS(MNUMS)
      REAL                PSETUP(MSETUP),           PSHAPE(MSHAPE)
      REAL                PSURF(MSURF),             PTRIAD(MTRIAD)
      REAL                PWCAP(MWCAP),             PWIND(MWIND)
      REAL                SIGMAG
      REAL                SPPARM(MSPPAR),SY0,       U10
      REAL                WDIC,          WDIP,      HSRERR
      REAL                PQUAD(MQUAD)
      REAL                RCOMPT(50000,5)
      REAL                PDIFFR(MDIFFR)                                  40.21
      LOGICAL             ACUPDA
      LOGICAL             BNDCHK,      BNAUT,       ONED,        BRESCL
      LOGICAL             OFFSRC                                          40.80
      LOGICAL             LADDS                                           40.85
      LOGICAL             CSETUP
!$OMP THREADPRIVATE(ICMAX,CSETUP)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE SWCOMM3

      MODULE SWCOMM4
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file SWCOMM4.INC
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     ---
!
!  7. Local variables
!
!     *** information for test output ***
!
! ICOTES [     0] minimum value for ITEST,
!                 =30; for command COTES, not documented in the manual
!                 =cotes; for command COTES [cotes], not documented
! INTES  [     0] testing parameter,
!                 =30; for command INTE, not documented in the manual
!                 =intes; for command INTE [intes], not documented
! IOUTES [     0] minimum value for ITEST,
!                 =30; for command OUTE, not documented in the manual
!                 =itest; for command OUTE [itest], not documented
! IPTST  [      ] sequence number of a test point
! IFPAR  [     0] unit ref. number for output of parameters in test points
! IFS1D  [     0] unit ref. number for output of 1D spectra of source terms
!                 if used, value is made non-zero by subr FOR
! IFS2D  [     0] unit ref. number for output of 2D spectra of source terms
!                 if used, value is made non-zero by subr FOR
! LXDMP  [    -1] grid counter for a test point in the x-direction,
!                 =ix; set by command TEST ... POI <[ix] [iy]>
! LYDMP  [     0] grid counter for a test point in the y-direction,
!                 =iy; set by command TEST ... POI <[ix] [iy]>
! MAXMES [   200] not used
!                 =maxmes; set by command SET [maxmes]
! NEGMES [     0] not used
! NPTST  [     0] number of test points; set by command TEST
! NPTSTA [     1] number of test points, equal to MAX(1,NPTST)
! TESTFL [   FAL] test output must/must not be made, mainly for testpoints
! UNDFLW [1.E-15] small number to prevent underflows
!
      INTEGER ICOTES,      INTES,       IOUTES
      INTEGER IPTST
      INTEGER IFPAR,       IFS1D,       IFS2D
      INTEGER LXDMP,       LYDMP,       MAXMES,      NEGMES
      INTEGER NPTST,       NPTSTA
      REAL    UNDFLW
      LOGICAL TESTFL
!$OMP THREADPRIVATE(IPTST,TESTFL)
!
!     *** higher order propagation and spherical coordinates ***
!
! COSLAT [    10] cos of latitude; =1 for Cartesian coordinates
! KREPTX [     0] if >0, the domain repeats itself in x-direction (primarily intended for
!                 propagation around the globe)
! KSPHER [     0] indicates whether spherical coordinates are used, and which projection method
!                 0=Cartesian coordinates, >0=spherical coordinates
! LENDEG [   1E5] length of a degree of the sphere
! PROJ_METHOD[ 0] projection method; 0=(quasi-)Cartesian,
!                 1=uniform Mercator (only spherical coordinates)
! PROPFL [CALCUL] indicates whether flux-limiting in spectral space is used
! PROPSC [CALCUL] indicates which numerical scheme is to be used for spatial propagation.
!                 1=first order (BSBT), 2=SORDUP, 3=3rd order (S&L)
! PROPSS [     2] indicates which numerical scheme is to be used in stationary computations
!                 1=first order (BSBT), 2=SORDUP
! PROPSN [     3] indicates which numerical scheme is to be used in nonstationary computations
!                 1=first order (BSBT), 3=3rd order (S&L)
! PROPSL [CALCUL] indicates which numerical scheme is used locally
! REARTH [   6E6] radius of the earth
! WAVAGE [    0.] indicates "wave age" parameter (used in counteracting garden-sprinkler effect
!                 in subroutine SANDL)
!
      INTEGER PROPSC,    KSPHER,    KREPTX
      INTEGER PROPSL
      INTEGER PROPSS,    PROPSN,    PROJ_METHOD
      INTEGER PROPFL
      REAL    WAVAGE,    REARTH,    LENDEG
      REAL    COSLAT(10)
!$OMP THREADPRIVATE(COSLAT,PROPSL)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE SWCOMM4

      MODULE TIMECOMM
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file TIMECOMM.INC
!
!  2. Purpose
!
!     Common variables used by the subroutines in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!TIMG!
!TIMG!     MXTIMR [  10] maximum number of simultaneous timings, dimension of
!TIMG!                   work-arrays LISTTM, TIMERS
!TIMG!     NSECTM [ 300] number of sections that may be defined for an application
!TIMG!
!TIMG      INTEGER NSECTM, MXTIMR
!TIMG      PARAMETER (NSECTM=300, MXTIMR=10)
!
!  7. Local variables
!
!     *** Time related variables for the computation ***
!
! TINIC     [    ] Start time and date of the computation (in seconds since
!                  the reference day (REFDAY))
!                  =tbegc; set by command COMP [tbegc] ...
! DT        [    ] Time step of the computation (in seconds)
!                  =deltc; set by command COMP ... [deltc] ...
!                  =deltc*60; set by command COMP ... [deltc] MI ...
!                  =deltc*60*60; set by command COMP ... [deltc] HR ...
!                  =deltc*60*60*24; set by command COMP ... [deltc] DA ...
! TFINC     [    ] End time and date of the computation (in seconds since
!                  the reference day (REFDAY))
! TIMCO     [    ] Time and date of the computation during the simulation (in
!                  seconds since the reference day (REFDAY))
!
      REAL TINIC, DT, TFINC, TIMCO
!
!     *** Time related variables for nested runs ***
!
! BEGBOU    [    ] Start time for the non-stationary boundary conditions
!                  (in seconds since the reference day (REFDAY)) in the case
!                  of nested runs. Read from the nest file
! IFACMX    [   1] set by command BOU (STAT) NE ... [ifacmx] ...
!                  (NOT documented)
! IFACMY    [   1] set by command BOU (STAT) NE ... [ifacmy]
!                  (NOT documented)
! TIMERB    [    ] Last time that non-stationary boundary conditions has been read
!                  (in seconds since the reference day (REFDAY))) in the case
!                  of nested runs. Read from the nest file
! TINTBO    [    ] Time step between non-stationary boundary conditions (in seconds)
!                  in the case of nested runs. Read from the nest file
!
      INTEGER  IFACMX, IFACMY
      REAL     BEGBOU, TIMERB, TINTBO
!TIMG!
!TIMG!     *** contains the cpu and wall-clock times ***
!TIMG!
!TIMG! DCUMTM    [ 600] cumulative time; columns 1,2: cpu-time, wall-clock time
!TIMG! LASTTM    [   1] last occupied position in LISTTM, 0 if all positions in
!TIMG!                  LISTTM are free
!TIMG! LISTTM    [  10] list of section numbers for all running/active timers
!TIMG!                  in array TIMERS. A value of -1 signals that the
!TIMG!                  corresponding timer is not running
!TIMG! NCUMTM    [ 300] for each section of the application the number of timings
!TIMG!                  that contributed to time in DCPUTM
!TIMG! TIMERS    [ 600] start-time of active timers; columns 1,2: cpu-time,
!TIMG!                  wall-clock time
!TIMG!
!TIMG      INTEGER, SAVE          :: NCUMTM(NSECTM), LISTTM(MXTIMR), LASTTM
!TIMG      DOUBLE PRECISION, SAVE :: DCUMTM(NSECTM,2), TIMERS(MXTIMR,2)
!TIMG!$OMP THREADPRIVATE(DCUMTM,TIMERS,NCUMTM,LISTTM,LASTTM)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
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
!     ---
!
! 13. Source text
!
      END MODULE TIMECOMM
      MODULE COUPLE

         real, dimension( : ), allocatable :: BLKND, BLKNDC, OURQT
         integer, dimension( : ), allocatable :: CROSS, BGRIDP

         real, dimension( :,: ), allocatable :: COMPDA
         real, dimension( :,:,: ), allocatable :: AC1
         real, dimension( :,:,:,: ), allocatable :: BSPECS

      END MODULE COUPLE
