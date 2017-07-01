!
!     SWAN - SERVICE ROUTINES
!
!  Contents of this file
!
!     READXY
!     REFIXY
!     INFRAM
!     DISTR
!     KSCIP1
!     AC2TST
!     CVCHEK                                                              30.60
!     CVMESH                                                              30.60
!     NEWTON                                                              30.60
!     EVALF                                                               30.60
!     SWOBST                                                              30.60
!     SWOBSTO                                                             40.86
!     TCROSS                                                              40.04
!     SWTRCF
!     SSHAPE                                                              40.00
!     SINTRP                                                              40.00
!     HSOBND                                                              32.01
!     CHGBAS                                                              40.00
!     GAMMA                                                               40.00
!     WRSPEC                                                              40.00
!TIMG!     SWTSTA                                                              40.23
!TIMG!     SWTSTO                                                              40.23
!TIMG!     SWPRTI                                                              40.23
!     TXPBLA                                                              40.23
!     INTSTR                                                              40.23
!     NUMSTR                                                              40.23
!     SWCOPI                                                              40.23
!     SWCOPR                                                              40.23
!
!  functions:
!  ----------
!  DEGCNV  (converts from cartesian convention to nautical and            32.01
!           vice versa)                                                   32.01
!  ANGRAD  (converts radians to degrees)                                  32.01
!  ANGDEG  (converts degrees to radians)                                  32.01
!
!  subroutines:
!  ------------
!  HSOBND  (Hs is calculated after a SWAN computation at all sides.       32.01
!           The calculated wave height from SWAN is then compared with    32.01
!           the wave heigth as provided by the user                       32.01
!
!***********************************************************************
!                                                                      *
      SUBROUTINE READXY (NAMX, NAMY, XX, YY, KONT, XSTA, YSTA)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE SWCOMM2                                                         40.41
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
!     40.22: John Cazes and Tim Campbell
!     40.13: Nico Booij
!     40.51: Marcel Zijlema
!
!  1. UPDATE
!
!       Nov. 1996               offset values are added to standard values
!                               because they will be subtracted later
!     40.13, Nov. 01: a valid value for YY is required if a valid value
!                     for XX has been given; ocpcomm1.inc reactivated
!     40.51, Feb. 05: correction to location points equal to offset values
!
!  2. PURPOSE
!
!       Read x and y, initialize offset values XOFFS and YOFFS
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NAMX, NAMY   inp char    names of the two coordinates as given in
!                                the user manual
!       XX, YY       out real    values of x and y taking into account offset
!       KONT         inp char    what to be done if values are missing
!                                see doc. of INDBLE (Ocean Pack doc.)
!       XSTA, YSTA   inp real    standard values of x and y
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       INDBLE (Ocean Pack)
      LOGICAL EQREAL
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
!       ----------------------------------------------------------------
!       Read x and y in double prec.
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      CHARACTER  NAMX *(*), NAMY *(*), KONT *(*)
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'READXY')
!
      CALL INDBLE (NAMX, XTMP, KONT, DBLE(XSTA)+DBLE(XOFFS))
      IF (CHGVAL) THEN                                                    40.13
!       a valid value was given for XX                                    40.13
        CALL INDBLE (NAMY, YTMP, 'REQ', DBLE(YSTA)+DBLE(YOFFS))           40.13
      ELSE                                                                40.13
        CALL INDBLE (NAMY, YTMP, KONT, DBLE(YSTA)+DBLE(YOFFS))
      ENDIF                                                               40.13
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      IF (.NOT.EQREAL(XOFFS,REAL(XTMP))) THEN                             40.51
         XX = REAL(XTMP-DBLE(XOFFS))
      ELSE IF (OPTG.EQ.3) THEN                                            40.51
         XX = 1.E-5                                                       40.51
      ELSE                                                                40.51
         XX = 0.                                                          40.51
      END IF                                                              40.51
      IF (.NOT.EQREAL(YOFFS,REAL(YTMP))) THEN                             40.51
         YY = REAL(YTMP-DBLE(YOFFS))
      ELSE IF (OPTG.EQ.3) THEN                                            40.51
         YY = 1.E-5                                                       40.51
      ELSE                                                                40.51
         YY = 0.                                                          40.51
      END IF                                                              40.51
!
      RETURN
! * end of subroutine READXY  *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE REFIXY (NDS, XX, YY, IERR)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM2                                                         40.41
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
!     40.22: John Cazes and Tim Campbell
!     40.51: M. Zijlema
!
!  1. UPDATE
!
!       first version: 10.18 (Sept 1994)
!
!  2. PURPOSE
!
!       initialize offset values XOFFS and YOFFS, and shift XX and YY
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NDS          in  int     file reference number
!       XX, YY       out real    values of x and y taking into account offset
!       IERR         out int     error indicator: IERR=0: no error, =-1: end-
!                                of-file, =-2: read error
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
      LOGICAL EQREAL
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
!       ----------------------------------------------------------------
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      REAL             XX, YY
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'REFIXY')
!
      READ (NDS, *, END=10, ERR=20) XTMP, YTMP
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      IF (.NOT.EQREAL(XOFFS,REAL(XTMP))) THEN                             40.51
         XX = REAL(XTMP-DBLE(XOFFS))
      ELSE IF (OPTG.EQ.3) THEN                                            40.51
         XX = 1.E-5                                                       40.51
      ELSE                                                                40.51
         XX = 0.                                                          40.51
      END IF                                                              40.51
      IF (.NOT.EQREAL(YOFFS,REAL(YTMP))) THEN                             40.51
         YY = REAL(YTMP-DBLE(YOFFS))
      ELSE IF (OPTG.EQ.3) THEN                                            40.51
         YY = 1.E-5                                                       40.51
      ELSE                                                                40.51
         YY = 0.                                                          40.51
      END IF                                                              40.51
!
      IERR = 0
      RETURN
!     end of file
  10  IERR = -1
      RETURN
!     read error
  20  IERR = -2
      RETURN
! * end of subroutine REFIXY  *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION  INFRAM (XQQ, YQQ)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
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
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!       Checking whether a point given in frame coordinates is located
!       in the plotting frame (INFRAM = .TRUE.) or not (INFRAM = .FALSE.)
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       XQQ     REAL   input    X-coordinate (output grid) of the point
!       YQQ     REAL   input    Y-coordinate (output grid) of the point
!
!  5. SUBROUTINES CALLING
!
!       SPLSIT, PLNAME
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
!       ----------------------------------------------------------------
!       Give INFRAM initial value true
!       IF XQQ < 0, XQQ > XQLEN, YQQ < 0 OR YQQ > YQLEN, THEN
!           INFRAM = false
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'INFRAM')
!
      INFRAM = .TRUE.
      IF (XQQ .LT.    0.) INFRAM = .FALSE.
      IF (XQQ .GT. XQLEN) INFRAM = .FALSE.
      IF (YQQ .LT.    0.) INFRAM = .FALSE.
      IF (YQQ .GT. YQLEN) INFRAM = .FALSE.
!
      RETURN
! * end of function INFRAM *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE DISTR (CDIR, DIR, COEF, SPCDIR)                          20.43
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
!
!  0. Authors
!
!     30.82: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!      0.1 , Jul. 87: Standard heading added
!      0.2 , Dec. 89: Value for energy outside of the sector changed
!                     from 0. to 1.E-6
!            Oct. 90: Value for energy outside the sector changed to 1.E-10
!                     logical BDIR introduced to take care for case where
!                     none of the values is positive
!     30.82, Oct. 98: Updated description of several variables
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!       Computation of the distribution of the wave energy over the
!       sectors, according to the given directional spread.
!
!  3. METHOD
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
!
      REAL    SPCDIR(MDC,6)                                               30.82
!
!       CDIR    REAL   output   array containing the coefficients of
!                               (energy) distribution
!       DIR     REAL   input    main wave direction, in radians
!       COEF    REAL   input    coefficient of the directional distri-
!                               bution (cos**COEF)
!
!  5. SUBROUTINES CALLING
!
!       REINVA (HISWA/SWREAD), STARTB and SWIND (both HISWA/COMPU)
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. Common blocks used
!
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       For every direction of the grid do
!           If the direction deviates less than PI/2 from the main wave
!            direction, then
!               Compute the coefficient cos**n
!           Else
!               Coefficient is 1.E-10
!       -----------------------------------------------------------------
!       If any of the directions deviated less than PI/2
!       Then Compute the total of the coefficients
!            For every direction of the grid do
!                Divide the fraction of the distribution by the total
!       -----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL BDIR
      REAL CDIR(*)
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT, 'DISTR')
!
      SOMC = 0.
      BDIR = .FALSE.
      DO 10 ID0 = 1, MDC
        TETA = SPCDIR(ID0,1)
        ACOS = COS(TETA-DIR)
        IF (ACOS .GT. 0.) THEN
          BDIR = .TRUE.
          CDIR(ID0) = MAX (ACOS**COEF, 1.E-10)
          SOMC = SOMC + CDIR(ID0)
        ELSE
          CDIR(ID0) = 1.E-10
        ENDIF
  10  CONTINUE
      IF (BDIR) THEN
        CNORM = 1./(SOMC*DDIR)
        DO 20 ID0 = 1, MDC
          CDIR(ID0) = CDIR(ID0) * CNORM
  20    CONTINUE
      ENDIF
!
!     ***** test *****
!      IF(TESTFL .AND. ITEST .GE. 200)
       IF( ITEST .GE. 200)
     &WRITE(PRINTF,6010) CNORM,(CDIR(JJ) , JJ=1,MDC)
 6010   FORMAT (' Test DISTR',F10.3/(10E12.4))
!
      RETURN
! * end of subroutine DISTR *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE KSCIP1 (MMT, SIG, D, K, CG, N, ND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
!
      IMPLICIT NONE
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
!     30.81: Annette Kieftenburg
!     40.22: John Cazes and Tim Campbell
!     40.41: Marcel Zijlema
!     41.16: Greg Wilson
!
!  1. Updates
!
!     Aug. 94, ver. 10.10: arguments N and ND added
!     Dec. 98, ND corrected, argument list adjusted and IMPLICIT NONE added
!     40.41, Aug. 04: tables replaced by Pade and other formulas
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.16, Mar. 11: correction: add dk/dh to dn/dh
!
!  2. Purpose
!
!     Calculation of the wave number, group velocity, group number N
!     and the derivative of N w.r.t. depth (=ND)
!
!  3. Method
!
!     --
!
!  4. Argument variables
!
!     MMT     input    number of frequency points
!
      INTEGER   MMT
!
!     CG      output   group velocity
!     D       input    local depth
!     K       output   wave number
!     N       output   ratio of group and phase velocity
!     ND      output   derivative of N with respect to D
!     SIG     input    rel. frequency for which wave parameters
!                      must be determined
!
      REAL      CG(MMT), D,
     &          K(MMT), N(MMT), ND(MMT), SIG(MMT)
!
!  6. Local variables
!
!     C         phase velocity
!     FAC1      auxiliary factor
!     FAC2      auxiliary factor
!     FAC3      auxiliary factor
!     IENT      number of entries
!     IS        counter in frequency (sigma-space)
!     KND       dimensionless wave number
!     ROOTDG    square root of D/GRAV
!     WGD       square root of GRAV*D
!     SND       dimensionless frequency
!     SND2      = SND*SND
!
      INTEGER   IENT, IS
      REAL      KND, ROOTDG, SND, WGD, SND2, C, FAC1, FAC2, FAC3
!
!  8. Subroutines used
!
!     --
!
!  9. Subroutines calling
!
!     SWOEXA, SWOEXF (Swan/Output)
!
! 10. Error messages
!
!     --
!
! 11. Remarks
!
!     --
!
! 12. Structure
!
!     -----------------------------------------------------------------
!      Compute non-dimensional frequency SND
!      IF SND >= 2.5, then
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND according to
!        deep water theory
!      ELSE IF SND =< 1.e-6
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND
!        according to extremely shallow water
!      ELSE
!        Compute wave number K, group velocity CGO and the ratio of
!        group and phase velocity N by Pade and other simple formulas.
!        Compute the derivative of N w.r.t. D = ND.
!     -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL STRACE (IENT, 'KSCIP1')
!
      ROOTDG = SQRT(D/GRAV)                                               30.81
      WGD    = ROOTDG*GRAV                                                30.81
      DO 200 IS = 1, MMT
!       SND is dimensionless frequency
        SND = SIG(IS) * ROOTDG
        IF (SND .GE. 2.5) THEN
!       ******* deep water *******
          K(IS)  = SIG(IS) * SIG(IS) / GRAV                               30.81
          CG(IS) = 0.5 * GRAV / SIG(IS)                                   30.81
          N(IS)  = 0.5
          ND(IS) = 0.
        ELSE IF (SND.LT.1.E-6) THEN
!       *** very shallow water ***                                        30.81
          K(IS)  = SND/D                                                  30.81
          CG(IS) = WGD
          N(IS)  = 1.
          ND(IS) = 0.
        ELSE
          SND2  = SND*SND                                                 40.41
          C     = SQRT(GRAV*D/(SND2+1./(1.+0.666*SND2+0.445*SND2**2       40.41
     &                                  -0.105*SND2**3+0.272*SND2**4)))   40.41
          K(IS) = SIG(IS)/C                                               40.41
          KND   = K(IS)*D                                                 40.41
          FAC1  = 2.*KND/SINH(2.*KND)                                     40.41
          N(IS) = 0.5*(1.+FAC1)                                           40.41
          CG(IS)= N(IS)*C                                                 40.41
          FAC2  = SND2/KND                                                40.41
          FAC3  = 2.*FAC2/(1.+FAC2*FAC2)                                  40.41
          FAC2  = -K(IS)*(2.*N(IS)-1.)/(2.*D*N(IS))                       41.16
          ND(IS)= FAC1*(0.5/D - K(IS)/FAC3 + FAC2*(0.5/K(IS) - D/FAC3))   41.16 40.41
        ENDIF
  200 CONTINUE
!
      RETURN
!     end of subroutine KSCIP1 *
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE AC2TST (XYTST, AC2,KGRPNT)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         41.13
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
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
!     0. Authors
!
!
!
!
!
      INTEGER   XYTST(*) ,KGRPNT(MXC,MYC)
      REAL      AC2(MDC,MSC,MCGRD)                                        30.21
!.................................................................
      IF ( ITEST .GE. 100 .AND. TESTFL) THEN
        DO II = 1, NPTST
          IF (OPTG.NE.5) THEN                                             40.80
             IX = XYTST(2*II-1)
             IY = XYTST(2*II)
             INDEX = KGRPNT(IX,IY)
             WRITE (PRINTF, 618) IX+MXF-2, IY+MYF-2, KGRPNT(IX,IY)        40.30
          ELSE                                                            40.80
             INDEX = XYTST(II)                                            40.80
             WRITE (PRINTF, 619) INDEX                                    40.80
          ENDIF                                                           40.80
          DO ID = 1, MDC
            WRITE (PRINTF, 620) (AC2(ID,IS,INDEX), IS=1,MIN(10,MSC))      30.21
          ENDDO
        ENDDO
      ENDIF
 618  FORMAT(/,'Spectrum for test point(index):', 2I5,2X,'(',I5,')')
 619  FORMAT(/,'Spectrum for test point: (',I5,')')                       40.80
 620  FORMAT (10(1X,E12.4))
      RETURN
      END
!****************************************************************
!
      SUBROUTINE CVCHEK (KGRPNT, XCGRID, YCGRID)                          30.72
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!            May  96: New subroutine
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.13, Mar. 01: messages corrected and extended
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Checks whether the given curvilinear grid is correct
!     also set the value of CVLEFT.
!
!  3. Method
!
!     Going around a mesh in the same direction the interior
!     of the mesh must be always in the same side if the
!     coordinates are correct
!
!  4. Argument variables
!
!     KGRPNT: input  Array of indirect addressing
!
      INTEGER KGRPNT(MXC,MYC)                                             30.72
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!
!     5. SUBROUTINES CALLING
!
!        SWRBC
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!     FIRST = True
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy),
!                    K3 = KGRPNT(ix+1,iy+1)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1),
!                    K3 = KGRPNT(ix,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1),
!                    K3 = KGRPNT(ix,iy)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy),
!                    K3 = KGRPNT(ix+1,iy)
!                 ---------------------------------------------------
!                 If K1>1 and K2>1 and K3>1
!                 Then Det = (xpg(K3)-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (ypg(K3)-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If FIRST
!                      Then Make FIRST = False
!                           If Det>0
!                           Then Make CVleft = False
!                           Else Make CVleft = True
!                      ----------------------------------------------
!                      If ((CVleft and Det<0) or (not CVleft and Det>0))
!                      Then Write error message with IX, IY, ISIDE
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!****************************************************************
!
!
      LOGICAL  FIRST
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'CVCHEK')
!
!     test output
!
      IF (ITEST .GE. 150 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,186)
 186    FORMAT(/,' ... Subroutine CVCHEK...',
     &  /,2X,'POINT( IX, IY),  INDEX,      COORDX,       COORDY')
        ICON = 0
        DO 5 IIY = 1, MYC
          DO 6 IIX = 1, MXC
            ICON = ICON + 1
            WRITE(PRINTF,7)IIX-1,IIY-1,KGRPNT(IIX,IIY),
     &      XCGRID(IIX,IIY)+XOFFS, YCGRID(IIX,IIY)+YOFFS                  30.72 40.13
 6        CONTINUE
 5      CONTINUE
      ENDIF
 7    FORMAT(4X,I5,1X,I5,3X,I4,5X,F10.2,4X,F10.2)
!
      FIRST = .TRUE.
!
      DO 10 IX = 1,MXC-1
        DO 15 IY = 1,MYC-1
          DO 20 ISIDE = 1,4
            IF (ISIDE .EQ. 1) THEN
              IX1 = IX                                                    40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 2) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 3) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY                                                    40.13
            ELSE IF (ISIDE .EQ. 4) THEN
              IX1 = IX                                                    40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY                                                    40.13
            ENDIF
            K1  = KGRPNT(IX1,IY1)                                         40.13
            XC1 = XCGRID(IX1,IY1)                                         40.13 30.72
            YC1 = YCGRID(IX1,IY1)                                         40.13 30.72
            K2  = KGRPNT(IX2,IY2)                                         40.13
            XC2 = XCGRID(IX2,IY2)                                         40.13 30.72
            YC2 = YCGRID(IX2,IY2)                                         40.13 30.72
            K3  = KGRPNT(IX3,IY3)                                         40.13
            XC3 = XCGRID(IX3,IY3)                                         30.72
            YC3 = YCGRID(IX3,IY3)                                         30.72
            DET   = 0.
            IF (K1 .GE. 2 .AND. K2 .GE. 2 .AND. K3 .GE. 2) THEN
              DET = ((XC3 - XC1) * (YC2 - YC1)) -
     &              ((YC3 - YC1) * (XC2 - XC1))
              IF (DET .EQ. 0.) THEN
!               three grid points on one line                             40.13
                CALL MSGERR (2,'3 comp. grid points on one line')         40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
 112            FORMAT (3(1X, 2I3, 2(1X, F14.4)))                         40.13
              ENDIF
!
              IF (FIRST) THEN
                FIRST = .FALSE.
                IF (DET .GT. 0.) THEN
                  CVLEFT = .FALSE.
                ELSE
                  CVLEFT = .TRUE.
                ENDIF
              ENDIF
              IF (     (      CVLEFT .AND. DET .GT. 0.)
     &            .OR. (.NOT. CVLEFT .AND. DET .LT. 0.)) THEN
!               crossing grid lines in a mesh                             40.13
                CALL MSGERR (2,'Grid angle <0 or >180 degrees')           40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
              ENDIF
            ENDIF
 20       CONTINUE
 15     CONTINUE
 10   CONTINUE
      RETURN
!     *** end of subroutine CVCHEK ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID, KGRBND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
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
!     30.72: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, May  98: procedure for points outside grid accelerated
!     40.00, Feb  99: procedure extended for 1D case
!                     XOFFS and YOFFS added in write statements
!     40.02, Mar. 00: Fixed bug that placed dry testpoints outside computational grid
!     40.13, Mar. 01: message "CVMESH 2nd attempt .." suppressed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Nov. 04: search for boundary points improved
!
!  2. Purpose
!
!     procedure to find location in curvilinear grid for a point
!     given in problem coordinates
!
!  3. Method
!
!     First attempt: use Newton-Raphson method to find XC and YC
!     (Note: in the program XC and YC indicate the mesh and position in
!     the mesh) in a few steps; this may be most efficient if a series of
!     points is processed, because the previous point provides a good
!     first estimate.
!     This procedure may fail if the number of iterations is larger than
!     a previously set limit (default=5).
!
!     If the first attempt fails then determine whether the points (XP,YP)
!     is inside the mesh. If so, then the Newton-Raphson procedure is used
!     again with the pivoting point like first guess. Otherwise, scan the
!     boundaries whether the point is on the boundaries. If this fails, it
!     may be concluded that the point (XP,YP) is outside the grid.
!
!  4. Argument variables
!
!     XCGRID  input  Coordinates of computational grid in x-direction     30.72
!     YCGRID  input  Coordinates of computational grid in y-direction     30.72
!     XP, YP  input  a point given in problem coordinates
!     XC, YC  outp   same point in computational grid coordinates
!
      REAL     XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                        30.72
      REAL     XP, YP, XC, YC
!
!     KGRPNT   input   array(MXC,MYC)  grid numbers
!                      if KGRPNT <= 1, point is not in comp. grid.
!     KGRBND   input   lists all boundary grid points consecutively
!
      INTEGER  KGRPNT(MXC,MYC), KGRBND(*)                                 40.00
!
!     Local variables
!
!     MXITNR   number of iterations in Newton-Raphson procedure
!     IX, IY   counter of computational grid point
!     K1       address of grid point
!     IXMIN    counter of grid point closest to (XP,YP)
!     IYMIN    counter of grid point closest to (XP,YP)
!     IBND     counter of boundary grid points
!
      INTEGER       :: IX, IY, K1, IXMIN, IYMIN, IBND
      INTEGER, SAVE :: MXITNR = 0
      INTEGER, SAVE :: IENT = 0
!
!     INMESH   if True, point (XP,YP) is inside the computational grid
!     FINDXY   if True, Newton-Raphson procedure succeeded
!     ONBND    if True, given point is on boundary                        40.41
!
      LOGICAL  INMESH ,FINDXY, ONBND
!
!     DISMIN   minimal distance found                                     40.41
!     XPC1     user coordinate of a computational grid point
!     YPC1     user coordinate of a computational grid point
!     XC0      grid coordinate of grid point closest to (XP,YP)
!     YC0      grid coordinate of grid point closest to (XP,YP)
!
      REAL       :: DISMIN                                                40.41
      REAL       :: XPC1, YPC1, XC0, YC0
!
!  5. SUBROUTINES CALLING
!
!     SINCMP
!
!  6. SUBROUTINES USED
!
!       NEWTON
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
!     --------------------------------------------------------------
!     Determine XC and YC from XP and YP using Newton-Raphson iteration
!     process
!     If (XC and YC were found) then
!       Procedure is ready; Return values of XC and YC
!       return
!     else
!     ---------------------------------------------------------------------
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             Inmesh = True
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy)
!                 ----------------------------------------------------------
!                 If K1>0 and K2>0
!                 Then Det = (xp-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (yp-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If ((CVleft and Det>0) or (not CVleft and Det<0))
!                      Then Make Inmesh = False
!                      Else  Inmesh = true and XC = IX and YC = IY
!                 Else Make Inmesh = False
!             --------------------------------------------------------
!             If Inmesh
!             Then Determine XC and YC using Newton-Raphson iteration
!                  process
!                  Procedure is ready; Return values of XC and YC
!     ---------------------------------------------------------------------
!     No mesh is found: Make XC and YC = exception value
!     Return values of XC and YC
!     ---------------------------------------------------------------------
!
!****************************************************************
!
!
      IF (LTRACE) CALL STRACE (IENT,'CVMESH')
!
      IF (ONED) THEN
        CALL NEWT1D  (XP, YP, XCGRID, YCGRID, KGRPNT,                     40.00
     &                XC ,YC ,FINDXY)
        IF (.NOT.FINDXY) THEN
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
          ENDIF
        ENDIF
        GOTO 99
      ELSE
!       two-dimensional computation
        XC = 1.
        YC = 1.
!       --- First attempt, to find XC,YC with Newton-Raphson method
        MXITNR = 5
        CALL NEWTON  (XP, YP, XCGRID, YCGRID,                             40.00
     &                MXITNR ,ITER, XC ,YC ,FINDXY)                       40.41 40.02
        IF ((ITEST .GE. 150 .OR. INTES .GE. 20) .AND. FINDXY) THEN        40.02
           WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC                    40.03
        ENDIF
 25     FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
        IF (FINDXY) GOTO 80
!
        IF (INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)) THEN
!         --- select grid point closest to (XP,YP)
          DISMIN = 1.E20
          DO 50 IX = 1,MXC
            DO 40 IY = 1,MYC
              K1  = KGRPNT(IX,IY)
              IF (K1.GT.1) THEN
                XPC1 = XCGRID(IX,IY)
                YPC1 = YCGRID(IX,IY)
                DISXY = SQRT ((XP-XPC1)**2 + (YP-YPC1)**2)
                IF (DISXY .LT. DISMIN) THEN
                  IXMIN  = IX
                  IYMIN  = IY
                  DISMIN = DISXY
                ENDIF
              ENDIF
  40        CONTINUE
  50      CONTINUE
!         second attempt using closest grid point as first guess
          MXITNR = 20
          XC0 = REAL(IXMIN)
          YC0 = REAL(IYMIN)
!         ITEST condition changed from 20 to 120                          40.13
          IF (ITEST.GE.120) WRITE (PRTEST, 55) XP+XOFFS ,YP+YOFFS ,       40.13
     &          XC0-1. ,YC0-1.
  55      FORMAT (' CVMESH 2nd attempt, (XP,YP)=','(',F12.4,',',F12.4,
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
          DO KORNER = 1, 4
            IF (KORNER.EQ.1) THEN
              XC = XC0 + 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.2) THEN
              XC = XC0 - 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.3) THEN
              XC = XC0 - 0.2
              YC = YC0 - 0.2
            ELSE
              XC = XC0 + 0.2
              YC = YC0 - 0.2
            ENDIF
            CALL NEWTON  (XP, YP, XCGRID, YCGRID,                         40.00
     &                    MXITNR ,ITER, XC ,YC ,FINDXY)                   40.41 40.02
            IF (FINDXY) THEN
              IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC               40.00
              ENDIF
              GOTO 80
            ENDIF
          ENDDO
          IF (ITER.GE.MXITNR) THEN                                        40.41
             WRITE (PRINTF, 75) XP+XOFFS, YP+YOFFS, MXITNR                40.00
  75         FORMAT (' search for point with location ', 2F12.4,          40.41
     &               ' fails in', I3, ' iterations')                      40.41
          END IF                                                          40.41
        ELSE
!         scan boundary to see whether the point is close to the boundary
          DISMIN=99999.                                                   40.41
          ONBND =.FALSE.                                                  40.41
          IX1 = 0                                                         40.41
          IY1 = 0                                                         40.51
          IX2 = 0
          DO IBND = 1, NGRBND
            IF (IX2.NE.0) THEN                                            40.41
               IX1 = IX2
               IY1 = IY2
               XP1 = XP2
               YP1 = YP2
            END IF
            IX2 = KGRBND(2*IBND-1)
            IY2 = KGRBND(2*IBND)
            IF (IX2.NE.0 .AND. (ABS(IX2-IX1).GT.1 .OR.                    40.51
     &                          ABS(IY2-IY1).GT.1)) IX1 = 0               40.51
            IF (IX2.GT.0) THEN
              XP2 = XCGRID(IX2,IY2)
              YP2 = YCGRID(IX2,IY2)
              IF (IBND.GT.1 .AND. IX1.GT.0) THEN                          40.51
!               --- determine relative distance from boundary segment
!                   with respect to the length of that segment
                SLEN2  = (XP2-XP1)**2 + (YP2-YP1)**2
                RELDIS = ABS((XP-XP1)*(YP2-YP1)-(YP-YP1)*(XP2-XP1)) /
     &                   SLEN2
                IF (RELDIS.LT.0.01) THEN                                  40.41
!                 --- determine location on the boundary section
                  IF (RELDIS-DISMIN.LE.0.01) THEN                         40.41
                     DISMIN = RELDIS                                      40.41
                     RELLOC = ((XP-XP1)*(XP2-XP1)+(YP-YP1)*(YP2-YP1)) /
     &                        SLEN2
                     IF (RELLOC.GE.-0.001 .AND. RELLOC.LE.1.001) THEN     40.41
                        RELLCM = RELLOC                                   40.41
                        IF (RELLCM.LT.0.01) RELLCM=0.                     40.41
                        IF (RELLCM.GT.0.99) RELLCM=1.                     40.41
                        IX1M  = IX1                                       40.41
                        IX2M  = IX2                                       40.41
                        IY1M  = IY1                                       40.41
                        IY2M  = IY2                                       40.41
                        ONBND = .TRUE.                                    40.41
                     ENDIF                                                40.41
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          IF (ONBND) THEN                                                 40.41
             XC = FLOAT(IX1M) + RELLCM * FLOAT(IX2M-IX1M) - 1.            40.41
             YC = FLOAT(IY1M) + RELLCM * FLOAT(IY2M-IY1M) - 1.            40.41
             IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF, 65) XP+XOFFS, YP+YOFFS, XC, YC              40.00
  65            FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &                  ') is on the boundary, (XC,YC)=(',
     &                  F9.2,',',F9.2,')')
             ENDIF
             GOTO 80
          ENDIF
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
  85        FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &              ') is outside grid')
          ENDIF
          GOTO 99
        ENDIF
      ENDIF                                                               40.00
  80  IF (KGRPNT(INT(XC+3.001)-2,INT(YC+3.001)-2).LE.1) THEN
         WRITE (PRINTF, 90) XP+XOFFS, YP+YOFFS                            40.41
  90     FORMAT (' point with location ',2F12.4,' is not active')         40.41
         XC = -99.
         YC = -99.
      ENDIF
  99  RETURN
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.41
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
!     Nico Booij
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!       New function for curvilinear version (ver. 40.00). May '98
!       40.03, Dec 99: test output added; commons swcomm2 and ocpcomm4 added
!       40.41, Oct. 04: common blocks replaced by modules, include files removed
!       40.41, Nov. 04: search for points restricted to subdomain
!       40.51, Feb. 05: determining number of crossing points improved
!
!  2. Purpose
!
!       procedure to find whether a given location is
!       in the (curvilinear) computational grid
!
!  3. Method  suggested by Gerbrant van Vledder
!
!       draw a line from the point (XP,YP) in vertical direction
!       determine the number of crossings with the boundary of the
!       grid; if this number is even the point is outside
!
!  4. Argument variables
!
!
!     KGRBND   int  input   array containing boundary grid points
!
      INTEGER  KGRBND(*)
!
!     XP, YP    real, input   a point given in problem coordinates
!     XCGRID    real, input   array(IX,IY) x-coordinate of a grid point
!     YCGRID    real, input   array(IX,IY) y-coordinate of a grid point
!
      REAL     XCGRID(MXC,MYC) ,YCGRID(MXC,MYC),
     &         XP, YP
!
!  5. Parameter variables
!
!  6. Local variables
!
!     NUMCRS   number of crossings with boundary outline
!
      INTEGER  NUMCRS, IX1, IY1, IX2, IY2
      REAL     XP1, XP2, YP1, YP2, YPS, RELDIS, RELDO
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!       CVMESH
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     --------------------------------------------------------------
!     numcros = 0
!     For all sections of the boundary do
!         determine coordinates of end points (XP1,YP1) and (XP2,YP2)
!         If (XP1<XP and XP2>XP) or (XP1>XP and XP2<XP)
!         then If not (YP1<YP and YP2<YP)
!                   if YPS>YP
!                   then numcros = numcros + 1
!     ---------------------------------------------------------------
!     If numcros is even
!     Then Inmesh = False
!     Else Inmesh = True
!     ---------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE (IENT,'INMESH')
!
      IF (XP.LT.XCLMIN .OR. XP.GT.XCLMAX .OR.                             40.41
     &    YP.LT.YCLMIN .OR. YP.GT.YCLMAX) THEN                            40.41
        IF (ITEST.GE.70) WRITE (PRTEST, 22) XP+XOFFS, YP+YOFFS,           40.03
     &    XCLMIN+XOFFS, XCLMAX+XOFFS, YCLMIN+YOFFS, YCLMAX+YOFFS          40.41
  22    FORMAT (1X, 2F12.4, ' is outside region ', 4F12.4)
        INMESH = .FALSE.
        GOTO 90
      ENDIF
!
      IF (NGRBND.LE.0) THEN
        CALL MSGERR (3, 'grid outline not yet determined')
        RETURN
      ENDIF
!
      NUMCRS = 0
      IX1    = 0
      IY1    = 0
      IX2    = 0
      RELDIS = -1.

!     loop over the boundary of the computational grid
      DO IBND = 1, NGRBND
        IF (IX2.NE.0) THEN
           IX1 = IX2
           IY1 = IY2
           XP1 = XP2
           YP1 = YP2
        END IF
        IX2 = KGRBND(2*IBND-1)
        IY2 = KGRBND(2*IBND)
        IF (IX2.NE.0 .AND. (ABS(IX2-IX1).GT.1 .OR.                        40.51
     &                      ABS(IY2-IY1).GT.1)) IX1 = 0                   40.51
        IF (IX2.GT.0) THEN
          XP2 = XCGRID(IX2,IY2)
          YP2 = YCGRID(IX2,IY2)
          IF (ITEST.GE.180) WRITE (PRTEST, 28) XP2+XOFFS,                 40.03
     &    YP2+YOFFS
  28      FORMAT (' boundary point ', 2F12.4)
          IF (IBND.GT.1 .AND. IX1.GT.0) THEN                              40.51
            IF (((XP1.GT.XP).AND.(XP2.LE.XP)).OR.
     &          ((XP1.LE.XP).AND.(XP2.GT.XP))) THEN
              IF (YP1.GT.YP .OR. YP2.GT.YP) THEN
!               determine y-coordinate of crossing point
                YPS = YP1 + (XP-XP1) * (YP2-YP1) / (XP2-XP1)
!               determine relative distance from boundary segment         40.51
!               with respect to the length of that segment                40.51
                RELDO  = RELDIS                                           40.51
                RELDIS = ABS(YP-YPS) / SQRT((XP2-XP1)**2 + (YP2-YP1)**2)  40.51
                IF (YPS.GT.YP.AND.ABS(RELDIS-RELDO).GT.0.1) THEN          40.51
                  NUMCRS = NUMCRS + 1
                  IF (ITEST.GE.70) WRITE (PRTEST, 32) NUMCRS,             40.03
     &            XP+XOFFS, YP+YOFFS, YPS+YOFFS                           40.03
  32              FORMAT (' crossing ', I1, ' point ', 3F12.4)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!     point is inside the grid is number of crossings is odd
      IF (MOD(NUMCRS,2) .EQ. 1) THEN
        INMESH = .TRUE.
      ELSE
        INMESH = .FALSE.
      ENDIF
  90  RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWTON (XP, YP, XCGRID, YCGRID,                          40.00
     &                   MXITNR, ITER, XC, YC, FIND)                      40.41 40.02
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     30.80, Oct. 98: computation of update of XC,YC modified to avoid
!                     division by 0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Solve eqs. and find a point  (XC,YC) in a curvilinear grid (compt.
!     grid) for a given point (XP ,YP) in a cartesian grid (problem coord).
!
!  3. Method
!
!     In this subroutine the next equations are solved :
!
!                  @XP             @XP
!     XP(xc,yc) +  --- * @XC   +   --- * @YC  - XP(x,y) = 0
!                  @XC             @YC
!
!                  @YP             @YP
!     YP(xc,yc) +  --- * @XC   +    --- * @YC  - YP(x,y) = 0
!                  @XC             @YC
!
!     In the subroutine, next notation is used for the previous eqs.
!     XVC       + DXDXC * DXC   + DXDYC * DYC - XP  = 0.
!     YVC       + DYDXC * DXC   + DYDYC * DYC - YP  = 0.
!
!
!  4. Argument variables
!
! i   MXITNR: Maximum number of iterations                                30.82
!
      INTEGER MXITNR, ITER                                                40.41 40.02
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!  6. SUBROUTINES USED
!
!     STRACE
!
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
!       -----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0
      IF (LTRACE) CALL STRACE (IENT,'NEWTON')
!
      DXC    = 1000.
      DYC    = 1000.
      TOLDC  = 0.001
      FIND   = .FALSE.
!
      IF (ITEST .GE. 200) THEN
        WRITE(PRINTF,*) ' Coordinates in subroutine NEWTON '
        DO J = 1, MYC
          DO I = 1, MXC
            WRITE(PRINTF,30) I ,J ,XCGRID(I,J) ,YCGRID(I,J)               30.72
          ENDDO
        ENDDO
      ENDIF
 30   FORMAT(2(2X,I5),2(2X,E12.4))
!
      DO 14 K = 1 ,MXITNR
        ITER = K                                                          40.41
        I1   = INT(XC)                                                    40.00
        J1   = INT(YC)
        IF (I1 .EQ. MXC) I1 = I1 - 1
        IF (J1 .EQ. MYC) J1 = J1 - 1
        I2  = I1 + 1
        J2  = J1 + 1
        FJ1 = FLOAT(J1)
        FI1 = FLOAT(I1)
        FJ2 = FLOAT(J2)
        FI2 = FLOAT(I2)
!
        XVC   = (YC-FJ1)*((XC-FI1)*XCGRID(I2,J2)  +
     &                    (FI2-XC)*XCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*XCGRID(I2,J1)  +
     &                    (FI2-XC)*XCGRID(I1,J1))
        YVC   = (YC-FJ1)*((XC-FI1)*YCGRID(I2,J2)  +
     &                    (FI2-XC)*YCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*YCGRID(I2,J1)  +
     &                    (FI2-XC)*YCGRID(I1,J1))
        DXDXC = (YC -FJ1)*(XCGRID(I2,J2) - XCGRID(I1,J2)) +
     &          (FJ2-YC )*(XCGRID(I2,J1) - XCGRID(I1,J1))
        DXDYC = (XC -FI1)*(XCGRID(I2,J2) - XCGRID(I2,J1)) +
     &          (FI2-XC )*(XCGRID(I1,J2) - XCGRID(I1,J1))
        DYDXC = (YC -FJ1)*(YCGRID(I2,J2) - YCGRID(I1,J2)) +
     &          (FJ2-YC )*(YCGRID(I2,J1) - YCGRID(I1,J1))
        DYDYC = (XC -FI1)*(YCGRID(I2,J2) - YCGRID(I2,J1)) +
     &          (FI2-XC )*(YCGRID(I1,J2) - YCGRID(I1,J1))
!
        IF (ITEST .GE. 150)
     &    WRITE(PRINTF,35) K, XC-1., YC-1., XP, YP, XVC, YVC              40.00
 35     FORMAT(' NEWTON  iter=', I2, ' (XC,YC)=', 2(1X,F10.2),/,          40.00
     &         ' (XP,YP)=', 2(1X,F10.2),
     &         '  X,Y(XC,YC) = ', 2(1X,F10.2))
        IF (ITEST .GE. 180) WRITE(PRINTF,36)
     &     XCGRID(I1,J1), XCGRID(I1,J2), XCGRID(I2,J1), XCGRID(I2,J2),
     &     YCGRID(I1,J1), YCGRID(I1,J2), YCGRID(I2,J1), YCGRID(I2,J2),
     &                     DXDXC, DXDYC, DYDXC, DYDYC                     40.00
 36     FORMAT(' NEWTON grid coord:', 8(1x, F10.0), /
     &         '        deriv=', 4(1X,F10.2))                             40.00
!
!       *** the derivated terms of the eqs. are evaluated and  ***
!       *** the eqs. are solved                                ***
        DDEN = DXDXC*DYDYC - DYDXC*DXDYC                                  30.80
        DXP  = XP - XVC                                                   30.80
        DYP  = YP - YVC                                                   30.80
        IF ( DDEN.NE.0. ) THEN
           DXC = ( DYDYC*DXP - DXDYC*DYP) / DDEN                          30.80
           DYC = (-DYDXC*DXP + DXDXC*DYP) / DDEN                          30.80
        ENDIF
!
        XC = XC + DXC
        YC = YC + DYC
!
!       *** If the guess point (XC,YC) is outside of compt. ***
!       *** grid, put that point in the closest boundary    ***
        IF (XC .LT. 1. ) XC = 1.
        IF (YC .LT. 1. ) YC = 1.
        IF (XC .GT. MXC) XC = FLOAT(MXC)
        IF (YC .GT. MYC) YC = FLOAT(MYC)
!
        IF (ITEST .GE. 120 .OR. INTES .GE. 50 .OR. IOUTES .GE. 50)
     &    WRITE(PRINTF,42) DXC, DYC, XC-1., YC-1.                         40.00
 42     FORMAT(' (DXC,DYC)=', 2(1X,F10.2), ' (XC,YC)=', 2(1X,F10.2))      40.00
!
!       *** If the accuracy is reached stop the iteration,  ***
        IF (ABS(DXC) .LE. TOLDC .AND. ABS(DYC) .LE. TOLDC) THEN
!
          FIND = .TRUE.
          XC = XC -1.
          YC = YC -1.
          RETURN
        ENDIF
!
 14   CONTINUE
      RETURN
!     *** end of subroutine NEWTON ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWT1D (XP, YP, XCGRID, YCGRID, KGRPNT,                  40.00
     &                   XC, YC, FIND)
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
!  0. Authors
!
!     40.00, 40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New (adaptation from subr NEWTON for 1D case)
!     40.13, Feb. 01: DX and DY renamed to DELX and DELY (DX and DY are
!                     common var.); error in expression for RS corrected
!                     PRINTF replaced by PRTEST in test output
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Finds broken coordinate XC for a given point XP in a rectilinear grid
!
!  3. Method
!
!     In this subroutine the step on the computational grid is selected
!     for which
!
!           (X-X1).(X2-X1)
!     0 <= --------------- <= 1
!          (X2-X1).(X2-X1)
!
!     where X, X1 and X2 are vectors; X corresponds to (Xp,Yp)
!     X1 and X2 are two neighbouring grid points
!
!  4. Argument variables
!
! i   KGRPNT: Grid adresses                                               40.00
!
      INTEGER KGRPNT(MXC,MYC)                                             40.00
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!     Local variables:

      REAL :: DELX, DELY   ! grid line                                    40.13
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
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'NEWT1D')
!
      IF (ITEST .GE. 120) THEN                                            40.13
        WRITE(PRTEST,*) ' Coordinates in subroutine NEWT1D '              40.13
        DO I = 1, MXC
          WRITE(PRTEST,30) I, XCGRID(I,1)+XOFFS ,YCGRID(I,1)+YOFFS        40.13
        ENDDO
      ENDIF
 30   FORMAT(2X,I5,2(2X,E12.4))
!
      FIND = .FALSE.
      DO 40 IX = 2 ,MXC
        IF (KGRPNT(IX-1,1).GT.1) THEN
          X1 = XCGRID(IX-1,1)
          Y1 = YCGRID(IX-1,1)
        ELSE
          GOTO 40
        ENDIF
        IF (KGRPNT(IX,1).GT.1) THEN
          X2 = XCGRID(IX,1)
          Y2 = YCGRID(IX,1)
        ELSE
          GOTO 40
        ENDIF
!       both ends of the step are valid grid points
!       now verify whether projection of (Xp,Yp) is within the step
        DELX = X2 - X1                                                    40.13
        DELY = Y2 - Y1                                                    40.13
        RS = ((XP - X1) * DELX + (YP - Y1) * DELY) /                      40.13
     &              (DELX * DELX + DELY * DELY)                           40.13
        IF (RS.GE.0. .AND. RS.LE.1.) THEN
          FIND = .TRUE.
          XC = REAL(IX-2) + RS                                            40.00
          YC = 0.
          GOTO 50
        ENDIF
  40  CONTINUE
  50  RETURN
!     *** end of subroutine NEWT1D ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE EVALF (XC ,YC ,XVC ,YVC ,XCGRID ,YCGRID)                 30.72
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
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Evaluate the coordinates (in problem coordinates) of point (XC,YC)
!     given in computational coordinates
!
!  3. Method
!
!     Bilinear interpolation
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!       XC, YC      real, outp    point in computational grid coordinates
!       XVC, YCV    real, OUTP    same point  but in problem coordinates
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
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'EVALF')
!
      I  = INT(XC)
      J  = INT(YC)
!
!     *** If the guess point (XC,YC) is in the boundary   ***
!     *** where I = MXC or/and J = MYC the interpolation  ***
!     *** is done in the mesh with pivoting point         ***
!     *** (MXC-1, J) or/and (I,MYC-1)                     ***
!
      IF (I .EQ. MXC) I = I - 1
      IF (J .EQ. MYC) J = J - 1
      T = XC - FLOAT(I)
      U = YC - FLOAT(J)
!     *** For x-coord. ***
      P1 = XCGRID(I,J)                                                    30.72
      P2 = XCGRID(I+1,J)                                                  30.72
      P3 = XCGRID(I+1,J+1)                                                30.72
      P4 = XCGRID(I,J+1)                                                  30.72
      XVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
!     *** For y-coord. ***
      P1 = YCGRID(I,J)                                                    30.72
      P2 = YCGRID(I+1,J)                                                  30.72
      P3 = YCGRID(I+1,J+1)                                                30.72
      P4 = YCGRID(I,J+1)                                                  30.72
      YVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
      RETURN
!     *** end of subroutine EVALF ***
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOBST (XCGRID, YCGRID, KGRPNT, CROSS)                   40.31 30.70
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_OBSTA                                                         40.31

      IMPLICIT NONE                                                       40.04
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
!     30.70
!     30.72  IJsbrand Haagsma
!     30.74  IJsbrand Haagsma
!     40.04  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.31  Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.70, Feb. 98: check if neighbouring point is a true grid point
!                     loop over grid points moved from calling routine into this
!                     argument list changed
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.04, Nov. 99: IMPLICIT NONE added, header updated
!                   : Removed include files that are not used
!     40.28, Feb. 02: Adjustments for extended REFLECTION option
!     40.31, Oct. 03: changes w.r.t. obstacles
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Obtains all the data required to find obstacles and
!     use subroutine TCROSS to find them
!
!  3. Method
!
!  4. Argument variables
!
!     CROSS   output Array which contains 0's if there is no
!                    obstacle crossing
!                    if an obstacle is crossing between the
!                    central point and its neighbour CROSS is equal
!                    to the number of the obstacle
!     KGRPNT  input  Indirect addressing for computational grid points
!     XCGRID  input  Coordinates of computational grid in x-direction     30.72
!     YCGRID  input  Coordinates of computational grid in y-direction     30.72
!
      INTEGER KGRPNT(MXC,MYC), CROSS(2,MCGRD)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)                            30.72
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ICC     index
!     ICGRD   index
!     IENT    number of entries of this subroutine
!     ILINK   indicates which link is analyzed: 1 -> neighbour in x
!                                               2 -> neighbour in y
!     IX      counter of gridpoints in x-direction
!     IY      counter of gridpoints in y-direction
!     JJ      counter for number of obstacles
!     JP      counter for number of corner points of obstacles
!     NUMCOR  number of corner points of obstacle
!     X1, Y1  user coordinates of one end of grid link
!     X2, Y2  user coordinates of other end of grid link
!     X3, Y3  user coordinates of one end of obstacle side
!     X4, Y4  user coordinates of other end of obstacle side
!
      INTEGER    ICC, ICGRD, IENT, ILINK, IX, IY, JJ, JP
      INTEGER    NUMCOR
      REAL       X1, X2, X3, X4, Y1, Y2, Y3, Y4
      LOGICAL    XONOBST                                                  40.04
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  8. Subroutines used
!
!     TCROSS                                                              40.04
!     STRACE
!
      LOGICAL    TCROSS                                                   40.04
!
!  9. Subroutines calling
!
!     SWPREP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!       ----------------------------------------------------------------
!       Read number of obstacles from array OBSTA
!       For every obstacle do
!           Read number of corners of the obstacle
!           For every corner of the obstacle do
!               For every grid point do
!                   call function TCROSS to search if there is crossing   40.04
!                   point                                                 40.04
!                   between the line of two points of the stencil and the
!                   line of the corners of the obstacle.
!                   If there is crossing point then
!                   then CROSS(link,kcgrd) = number of the crossing obstacle
!                   else CROSS(link,kcgrd) = 0
!       ----------------------------------------------------------------
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWOBST')
!
      IF (NUMOBS .GT. 0) THEN
!       NUMOBS is the number of obstacles ***
        COBST => FOBSTAC                                                  40.31
        DO 120 JJ = 1, NUMOBS
!         number of corner points of the obstacle
          NUMCOR = COBST%NCRPTS
          IF (ITEST.GE. 120) THEN
            WRITE(PRINTF,50) JJ, NUMCOR
 50         FORMAT( ' Obstacle number : ', I4,'  has ',I4,' corners')
          ENDIF
!         *** X1 X2 X3 ETC. are the coordinates of point according ***
!         *** with the scheme in the subroutine TCROSS header      ***    40.04
          X3 = COBST%XCRP(1)                                              40.31
          Y3 = COBST%YCRP(1)                                              40.31
          IF (ITEST.GE. 120)  WRITE(PRINTF,30) 1,X3,Y3
          DO 110 JP = 2, NUMCOR
            X4 = COBST%XCRP(JP)                                           40.31
            Y4 = COBST%YCRP(JP)                                           40.31
            IF (ITEST.GE. 120) WRITE(PRINTF,30) JP,X4,Y4
  30        FORMAT(' Corner number:', I4,'    XP: ',E10.4,' YP: ',E11.4)
            DO 100 IX = 1, MXC
              DO 90 IY = 1, MYC
                ICC = KGRPNT(IX,IY)
                IF (ICC .GT. 1) THEN
                  X1 = XCGRID(IX,IY)                                      30.72
                  Y1 = YCGRID(IX,IY)                                      30.72
!
!                 *** "ILINK" indicates which link is analyzed. Initial  ***
!                 *** neighbour in x , second link with neighbouring in y***
                  DO 80 ILINK = 1, 2
                    IF (ILINK.EQ.1 .AND. IX.GT.1) THEN
                      X2    = XCGRID(IX-1,IY)                             30.72
                      Y2    = YCGRID(IX-1,IY)                             30.72
                      ICGRD = KGRPNT(IX-1,IY)                             30.70
                    ELSE IF (ILINK.EQ.2 .AND. IY.GT.1) THEN
                      X2    = XCGRID(IX,IY-1)                             30.72
                      Y2    = YCGRID(IX,IY-1)                             30.72
                      ICGRD = KGRPNT(IX,IY-1)                             30.70
                    ELSE
                      ICGRD = 0
                    ENDIF
                    IF (ICGRD.GT.1) THEN                                  30.70
!
!                     *** All links are analyzed in each point otherwise the   ***
!                     *** boundaries can be excluded                           ***
!
                      IF (TCROSS(X1, X2, X3, X4, Y1, Y2, Y3, Y4,          40.04
     &                           XONOBST)) THEN                           40.04
                        CROSS(ILINK,ICC) = JJ
                      ENDIF
                    ENDIF
  80              CONTINUE
                ENDIF
  90          CONTINUE
 100        CONTINUE
            X3 = X4
            Y3 = Y4
 110      CONTINUE
          IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT
          COBST => COBST%NEXTOBST
 120    CONTINUE
      ENDIF
!
      RETURN
! * end of subroutine SWOBST *
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWOBSTO (XCGRID, YCGRID, XP, YP, XC, YC, KGRPNT, CROSS,
     &                    MIP)
!                                                                      *
!***********************************************************************

      USE OCPCOMM4
      USE SWCOMM3
      USE SWCOMM4
      USE M_OBSTA

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Nico Booij                                    |
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
!     40.86: Nico Booij
!
!  1. Updates
!
!     40.86, Feb. 08: new subroutine, based on SWOBST
!
!  2. Purpose
!
!     Find out whether a line segment defined in user coordinates
!     crosses one or more obstacles
!     use subroutine TCROSS to find crossings
!
!  3. Method
!
!  4. Argument variables
!
!     CROSS   output false if no obstacle crossing
!                    true if an obstacle is crossing between the
!                    output point and neighbouring grid point
!     KGRPNT  input  Indirect addressing for computational grid points
!     MIP     input  number of output points
!     XC, YC  input  array containing indices of output points
!     XP, YP  input  user coordinates of output point
!     XCGRID  input  Coordinates of computational grid in x-direction
!     YCGRID  input  Coordinates of computational grid in y-direction
!
      INTEGER MIP
      INTEGER KGRPNT(MXC,MYC)
      LOGICAL CROSS(4,MIP)
      REAL    XC(MIP), YC(MIP), XP(MIP), YP(MIP)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)
!
!  5. Parameter variables
!
!  6. Local variables
!
      INTEGER :: IP   ! counter for number of output points
      INTEGER :: JJ   ! counter for number of obstacles
      INTEGER :: JP   ! counter for number of corner points of obstacles
      INTEGER :: NUMCOR  ! number of corner points of obstacle
      REAL :: X1, Y1  ! user coordinates of one end of line segment
      REAL :: X2, Y2  ! user coordinates of other end of line segment
      REAL :: X3, Y3  ! user coordinates of one end of obstacle side
      REAL :: X4, Y4  ! user coordinates of other end of obstacle side
      INTEGER :: JX(1:4), JY(1:4) ! grid counters for the 4 corners
      INTEGER :: JC            ! corner counter
      INTEGER :: JX1, JY1, JX2, JY2
      LOGICAL :: OUTSID
!
      TYPE(OBSTDAT), POINTER :: COBST
!
!  8. Subroutines used
!
!     STRACE
!
      LOGICAL :: TCROSS   ! determines whether two line segments cross
!
!  9. Subroutines calling
!
!     SWOUTP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
! 13. Source text
! ======================================================================

      LOGICAL :: XONOBST        ! not used
      INTEGER, SAVE :: IENT=0   ! number of entries of this subroutine
      IF (LTRACE) CALL STRACE (IENT,'SWOBSTO')
!
      CROSS = .FALSE.
!
      IF (NUMOBS .GT. 0) THEN
!       NUMOBS is the number of obstacles ***
        COBST => FOBSTAC
        DO JJ = 1, NUMOBS
!         number of corner points of the obstacle
          NUMCOR = COBST%NCRPTS
          IF (ITEST.GE. 120) THEN
            WRITE(PRINTF,50) JJ, NUMCOR
 50         FORMAT( ' Obstacle number : ', I4,'  has ',I4,' corners')
          ENDIF
!         *** X1 X2 X3 ETC. are the coordinates of point according ***
!         *** with the scheme in the subroutine TCROSS header      ***
          X3 = COBST%XCRP(1)
          Y3 = COBST%YCRP(1)
          IF (ITEST.GE. 120)  WRITE(PRINTF,30) 1,X3,Y3
          DO JP = 2, NUMCOR
             X4 = COBST%XCRP(JP)
             Y4 = COBST%YCRP(JP)
             IF (ITEST.GE. 120) WRITE(PRINTF,30) JP,X4,Y4
  30         FORMAT(' Corner number:', I4,'    XP: ',E10.4,
     &                                       ' YP: ',E11.4)
             DO IP = 1, MIP
                X1 = XP(IP)
                Y1 = YP(IP)
!
                IF (XC(IP) .LE. -0.5 .OR. YC(IP) .LE. -0.5) CYCLE
                OUTSID = .FALSE.
                JX1 = INT(XC(IP)+3.001) - 2
                JX2 = JX1 + 1
                IF (JX1.LT.0) OUTSID = .TRUE.
                IF (KREPTX .EQ. 0) THEN
                   IF (JX1.GT.MXC) OUTSID = .TRUE.
                   IF (JX1.EQ.MXC) JX2 = MXC
                   IF (JX1.EQ.0)   JX1 = 1
                ELSE
                   JX1 = 1 + MODULO (JX1-1,MXC)
                   JX2 = 1 + MODULO (JX2-1,MXC)
                ENDIF
                IF (ONED) THEN
                   JY1 = 1
                   JY2 = 1
                ELSE
                   JY1 = INT(YC(IP)+3.001) - 2
                   JY2 = JY1 + 1
                   IF (JY1.LT.0)   OUTSID = .TRUE.
                   IF (JY1.GT.MYC) OUTSID = .TRUE.
                   IF (JY1.EQ.MYC) JY2 = MYC
                   IF (JY1.EQ.0)   JY1 = 1
                ENDIF
                IF (.NOT.OUTSID) THEN
                   JX(1) = JX1
                   JY(1) = JY1
                   JX(2) = JX2
                   JY(2) = JY1
                   JX(3) = JX1
                   JY(3) = JY2
                   JX(4) = JX2
                   JY(4) = JY2
!
                   DO JC = 1, 4
                      IF (KGRPNT(JX(JC),JY(JC)).GT.1) THEN
                         X2 = XCGRID(JX(JC),JY(JC))
                         Y2 = YCGRID(JX(JC),JY(JC))
                         IF (TCROSS(X1, X2, X3, X4, Y1, Y2, Y3, Y4,
     &                              XONOBST)) CROSS(JC,IP) = .TRUE.
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
!
             X3 = X4
             Y3 = Y4
          ENDDO
          IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT
          COBST => COBST%NEXTOBST
        ENDDO
      ENDIF
!
      RETURN
!     end of subroutine SWOBSTO
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TCROSS (X1, X2, X3, X4, Y1, Y2, Y3, Y4, X1ONOBST)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE                                                       40.04
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
!     40.00  Gerbrant van Vledder
!     40.04  Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!       30.70, Feb 98: argument list simplified
!                      subroutine changed into logical function
!       40.00, Aug 98: division by zero prevented
!       40.04, Aug 99: method corrected, IMPLICIT NONE added, XCONOBST added,
!                      introduced TINY and EPSILON (instead of comparing to 0)
!                      replaced 0 < LMBD,MIU by  0 <= LMBD,MIU
!                      XCONOBST added to argument list
!       40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!       Find if there is an obstacle crossing the stencil in used
!
!  3. Method
!
!     For the next situation (A, B and C are the points in the stencil,
!     D and E  are corners of the obstacle
!
!
!      obstacle --> D(X3,Y3)
!                    *
!                     *
!                      *
!        (X2,Y2)        * (XC,YC)
!            B-----------@--------------------------A (X1,Y1)
!                        ^*                         /
!                   _____| *                       /
!                  |        *                     /
!                  |         *                   /
!         crossing point      *                 /
!                              *               /
!                               E             /
!                              (X4,Y4)       /
!                                           C
!
!
!       The crossing point (@) should be found solving the next eqs.
!       for LMBD and MIU.
!
!       | XC |    | X1 |           | X2 - X1 |
!       |    | =  |    | +  LMBD * |         |
!       | YC |    | Y1 |           | Y2 - Y1 |
!
!
!       | XC |    | X3 |           | X4 - X3 |                            40.04
!       |    | =  |    | +  MIU  * |         |
!       | YC |    | Y3 |           | Y4 - Y3 |                            40.04
!
!
!     If solution exist and (0 <= LMBD <= 1 and 0 <= MIU <= 1)            40.04
!     there is an obstacle crossing the stencil
!
!  4. Argument variables
!
!     X1, Y1  inp    user coordinates of one end of grid link
!     X2, Y2  inp    user coordinates of other end of grid link
!     X3, Y3  inp    user coordinates of one end of obstacle side
!     X4, Y4  inp    user coordinates of other end of obstacle side
!     X1ONOBST outp   boolean which tells whether (X1,Y1) is on obstacle
!
      REAL       EPS, X1, X2, X3, X4, Y1, Y2, Y3, Y4
      LOGICAL    X1ONOBST
!
!  5. Parameter variables
!
!  6. Local variables
!
!     A,B,C,D    dummy variables
!     DIV1       denominator of value of LMBD (or MIU)
!     E,F        dummy variables
!     IENT       number of entries
!     LMBD       coefficient in vector equation for stencil points (or obstacle)
!     MIU        coefficient in vector equation for obstacle (or stencil points)
!
      INTEGER    IENT
      REAL       A, B, C, D, DIV1, E, F, LMBD, MIU
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOBST
!     SWTRCF
!     OBSTMOVE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     Calculate MIU and LMBD
!     If 0 <= MIU, LMBD <= 1                                              40.04
!     Then TCROSS is .True.
!     Else TCROSS is .False.
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'TCROSS')
!
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))             40.04
      IF (EPS ==0.) EPS = TINY(X1)                                        40.04
      A    = X2 - X1
!     A not equal to zero
      IF (ABS(A) .GT. TINY(X1)) THEN                                      40.04
        B    = X4 - X3
        C    = X3 - X1
        D    = Y2 - Y1
        E    = Y4 - Y3
        F    = Y3 - Y1
      ELSE
!       exchange MIU and LMBD                                             40.04
        A    = X4 - X3
        B    = X2 - X1
        C    = X1 - X3
        D    = Y4 - Y3
        E    = Y2 - Y1
        F    = Y1 - Y3
      ENDIF
      DIV1 = ((A*E) - (D*B))                                              40.00
!
!     DIV1 = 0 means that obstacle is parallel to line through            40.04
!     stencil points, or (X3,Y3) = (X4,Y4);                               40.04
!     A = 0 means trivial set of equations X4= X3 and X2 =X1              40.04
!
      IF ((ABS(DIV1).LE.TINY(X1)) .OR.                                    40.04
     &    (ABS(A).LE.TINY(X1))) THEN                                      40.04
        MIU = -1.                                                         40.00
        LMBD = -1.                                                        40.04
      ELSE                                                                40.00
        MIU  = ((D*C) - (A*F)) / DIV1
        LMBD = (C + (B*MIU)) / A
      END IF                                                              40.00
!
      IF (MIU  .GE. 0. .AND. MIU  .LE. 1. .AND.                           40.04
     &    LMBD .GE. 0. .AND. LMBD .LE. 1.) THEN                           40.04
!
!       Only (X1,Y1) is checked, because of otherwise possible double     40.04
!       counting                                                          40.04
        IF ((LMBD.LE.EPS .AND. ABS(X2-X1).GT.EPS).OR.                     40.04
     &      (MIU .LE.EPS .AND. ABS(X2-X1).LE.EPS))THEN                    40.04
          X1ONOBST = .TRUE.                                               40.04
        ELSE                                                              40.04
          X1ONOBST = .FALSE.                                              40.04
        ENDIF                                                             40.04
!
!       *** test output ***
        IF (ITEST .GE. 120) THEN
          WRITE(PRINTF,70)X1,Y1,X2,Y2,X3,Y3,X4,Y4
  70      FORMAT(' Obstacle crossing  :',/,
     &    ' Coordinates of comp grid points and corners of obstacle:',/,
     &    ' P1(',E10.4,',',E10.4,')',' P2(',E10.4,',',E10.4,')',/,
     &    ' P3(',E10.4,',',E10.4,')',' P4(',E10.4,',',E10.4,')')
        ENDIF
!
        TCROSS = .TRUE.
      ELSE
        TCROSS = .FALSE.
      ENDIF
!
!     End of subroutine TCROSS
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      RECURSIVE SUBROUTINE OBSTMOVE (XCGRID, YCGRID, KGRPNT)              40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.81
      USE SWCOMM3                                                         40.41
      USE M_OBSTA                                                         40.31

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
!     40.09  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.31  Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.09, July 00: new subroutine
!     40.28, Feb. 02: Adjustments for extended REFLECTION option
!     40.31, Oct. 03: changes w.r.t. obstacles
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Move OBSTACLE points (X3,Y3) and (X4,Y4) a bit if computational gridcell
!     (X1,Y1) is on the OBSTACLE line piece.
!
!  3. Method
!
!     Add EPS*(dY,-dX) to OBSTACLE line piece coordinates so that movement of
!     these OBSTACLE points is perpendicular to the direction
!
!  4. Argument variables
!
!     KGRPNT  input  Indirect addressing for computational grid points
!     XCGRID  input  Coordinates of computational grid in x-direction
!     YCGRID  input  Coordinates of computational grid in y-direction
!
      INTEGER KGRPNT(MXC,MYC)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)
!
!  5. Parameter variables
!
!  6. Local variables
!
!     DISTA    distance between (X1,Y1) and (X2A,Y2A)
!     DISTB    distance between (X1,Y1) and (X2 ,Y2 )
!     DXA, DYA difference X1 - X2A respectively Y1 - Y2A
!     DXB, DYB difference X2 - X1 respectively Y2 - Y1
!     DXO, DYO difference X4 - X3 respectively Y4 - Y3
!     DXYO     distance between (X3,Y3) and (X4,Y4)
!     EPS      multiplication factor
!     ICC      index
!     ICGRD    index
!     IENT     number of entries of this subroutine
!     ILINK    indicates which link is analyzed: 1 -> neighbour in x
!                                                2 -> neighbour in y
!     IX       counter of gridpoints in x-direction
!     IY       counter of gridpoints in y-direction
!     JJ       counter for number of obstacles
!     JP       counter for number of corner points of obstacles
!     MOVED    boolean which tells whether OBSTACLE has been moved
!     NUMCOR   number of corner points of obstacle
!     X1, Y1   computational grid coordinates of one end of grid link
!     X2, Y2   computational grid coordinates of other end of grid link
!              i.e. neighbouring point of X1,Y1 associated with linknumber
!     X2A,Y2A  other neighbouring point of X1,Y1 associated with other
!              linknumber, or if invalid: = X2,Y2
!     X3, Y3   user coordinates of one end of obstacle side
!     X4, Y4   user coordinates of other end of obstacle side
!     XCONOBST boolean variable which tells whether XC in on OBSTACLE
!              line piece
!     XYEPS    displacement factor relative to local computational grid
!
      INTEGER    ICC, ICGRD, IENT, ILINK, IX, IY, JJ, JP
      INTEGER    NUMCOR
      REAL       DISTA, DISTB, DXA, DXB, DYA, DYB,
     &           DXO, DXYO, DYO, EPS, X1, X2, X2A, X3, X4,
     &           XYEPS, Y1, Y2, Y2A, Y3, Y4
      LOGICAL    MOVED, XCONOBST
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  8. Subroutines used
!
!     STRACE
!     TCROSS
!     MSGERR
!
      LOGICAL    TCROSS
!
!  9. Subroutines calling
!
!     SWPREP
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!       ----------------------------------------------------------------
!       Read number of obstacles from array OBSTA
!       For every obstacle do
!           Read number of corners of the obstacle
!           For every corner of the obstacle do
!               For every grid point do
!                   call function TCROSS to search if there is crossing point
!                   between the line of two points of the stencil and the
!                   line of the corners of the obstacle.
!                   If there is crossing point then
!                     If there computational gridpoint is on the obstacle
!                     then move corner points of obstacle perpendicular
!                     to obstacle with factor
!                     (XYEPS*DYO/DXYO,-XYEPS*DXO/DXYO)
!                     Moved = .True.
!                   If Moved then call OBSTMOVE again to check whether
!                   there are still computational grid points on OBSTACLE
!       ----------------------------------------------------------------
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'OBSTMOVE')
      MOVED = .FALSE.
      EPS =1.E-2
!
      IF (NUMOBS .GT. 0) THEN
        COBST => FOBSTAC                                                  40.31
        DO JJ = 1, NUMOBS
!         number of corner points of the obstacle
          NUMCOR = COBST%NCRPTS
          IF (ITEST.GE. 120) THEN
            WRITE(PRINTF,50) JJ, NUMCOR
 50         FORMAT( ' Obstacle number : ', I4,'  has ',I4,' corners')
          ENDIF
!         *** X1 X2 X3 ETC. are the coordinates of point according ***
!         *** with the scheme in the subroutine TCROSS header      ***
          X3 = COBST%XCRP(1)                                              40.31
          Y3 = COBST%YCRP(1)                                              40.31
          DO JP = 2, NUMCOR
            X4 = COBST%XCRP(JP)                                           40.31
            Y4 = COBST%YCRP(JP)                                           40.31
            IF (ITEST.GE. 120) WRITE(PRINTF,30) JP,X4,Y4
  30        FORMAT(' Corner number:', I4,'    XP: ',E10.4,' YP: ',E11.4)
            DO IX = 1, MXC
              DO IY = 1, MYC
                ICC = KGRPNT(IX,IY)
                IF (ICC .GT. 1) THEN
                  X1 = XCGRID(IX,IY)
                  Y1 = YCGRID(IX,IY)
!
!                 *** "ILINK" indicates which link is analyzed. Initial  ***
!                 *** neighbour in x , second link with neighbouring in y***
                  DO ILINK = 1, 2
                    IF (ILINK.EQ.1 .AND. IX.GT.1) THEN
                      X2    = XCGRID(IX-1,IY)
                      Y2    = YCGRID(IX-1,IY)
                      IF (IY.GT.1) THEN
                        X2A    = XCGRID(IX,IY-1)
                        Y2A    = YCGRID(IX,IY-1)
                      ELSE
                        X2A    = X2
                        Y2A    = Y2
                      ENDIF
                      ICGRD = KGRPNT(IX-1,IY)
                    ELSE IF (ILINK.EQ.2 .AND. IY.GT.1) THEN
                      X2    = XCGRID(IX,IY-1)
                      Y2    = YCGRID(IX,IY-1)
                      IF (IX.GT.1) THEN
                       X2A    = XCGRID(IX-1,IY)
                       Y2A    = YCGRID(IX-1,IY)
                      ELSE
                       X2A    = X2
                       Y2A    = Y2
                      ENDIF
                      ICGRD = KGRPNT(IX,IY-1)
                    ELSE
                     ICGRD = 0
                    ENDIF
                    IF (ICGRD.GT.1) THEN
!
!                     *** All links are analyzed in each point otherwise the   ***
!                     *** boundaries can be excluded                           ***
!
                      IF (TCROSS(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XCONOBST))THEN
                        IF (XCONOBST) THEN
                        DXA=(X1-X2A)
                        DYA=(Y1-Y2A)
                        DXB=(X2-X1)
                        DYB=(Y2-Y1)
                        DISTA=SQRT(DXA*DXA+DYA*DYA)
                        DISTB=SQRT(DXB*DXB+DYB*DYB)
                        XYEPS = EPS*MIN(DISTA,DISTB)
                        DXO = X4-X3
                        DYO = Y4-Y3
                        DXYO = SQRT(DXO*DXO+DYO*DYO)
!                       -DXO/DXYO and DYO/DXYO are used (instead of -DXO
!                       and DYO) because otherwise displacement is dependent
!                       on length of OBSTACLE line piece
                        COBST%XCRP(JP-1) = X3 + XYEPS*DYO/DXYO            40.31
                        COBST%YCRP(JP-1) = Y3 - XYEPS*DXO/DXYO            40.31
                        COBST%XCRP(JP  ) = X4 + XYEPS*DYO/DXYO            40.31
                        COBST%YCRP(JP  ) = Y4 - XYEPS*DXO/DXYO            40.31
                        CALL MSGERR (1, 'Obstacle points moved')
                        WRITE(PRINTF, 17) X3+XOFFS, Y3+YOFFS,             40.81
     &                        X4+XOFFS, Y4+YOFFS,                         40.81
     &                        X3+XYEPS*DYO/DXYO+XOFFS,                    40.81
     &                        Y3-XYEPS*DXO/DXYO+YOFFS,                    40.81
     &                        X4+XYEPS*DYO/DXYO+XOFFS,                    40.81
     &                        Y4-XYEPS*DXO/DXYO+YOFFS,                    40.81
     &                        X1+XOFFS,Y1+YOFFS                           40.81
  17                    FORMAT ('OBSTACLE POINTS (', F11.2, ',',  F11.2,
     &                         '), and (', F11.2,',',  F11.2,'),',
     &                         'moved to: (',  F11.2,',',
     &                         F11.2,'), and (', F11.2,',', F11.2,
     &                         '), because OBSTACLE line piece ',
     &                         'was on computational grid point (',
     &                         F11.2,',', F11.2,').')
                        X3 = X3 + XYEPS * DYO/DXYO
                        Y3 = Y3 - XYEPS * DXO/DXYO
                        X4 = X4 + XYEPS * DYO/DXYO
                        Y4 = Y4 - XYEPS * DXO/DXYO
                        MOVED = .TRUE.
                        ENDIF
                      ENDIF
                    ENDIF
                  END DO
                ENDIF
              END DO
            END DO
            X3 = X4
            Y3 = Y4
          END DO
          IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT
          COBST => COBST%NEXTOBST
        END DO
      ENDIF
      IF (MOVED) CALL OBSTMOVE(XCGRID, YCGRID, KGRPNT)                    40.31
      RETURN
! * end of subroutine OBSTMOVE *
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWTRCF (WLEV2 , CHS   ,                                  40.31 40.00
     &                   LINK  , OBREDF,                                  40.03
     &                   AC2   , REFLSO, KGRPNT, XCGRID,                  40.41 40.09
     &                   YCGRID, CAX,    CAY,    RDX,    RDY,    ANYBIN,  40.09
     &                   SPCSIG, SPCDIR)                                  40.13 40.28
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        40.41
      USE SWCOMM2                                                         40.66
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_OBSTA                                                         40.31
      USE M_PARALL                                                        40.31
      USE SwanGriddata                                                    40.80
      USE SwanCompdata                                                    40.80

      IMPLICIT NONE                                                       40.09
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
!     30.70
!     40.03  Nico Booij
!     40.08  Erick Rogers
!     40.09  Annette Kieftenburg
!     40.13  Nico Booij
!     40.14  Annette Kieftenburg
!     40.18  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.30  Marcel Zijlema
!     40.31  Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.66: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.70, Feb. 98: water level (WLEV2) replaced depth
!                     incident wave height introduced using argument
!                     CHS (sign. wave height in whole comput. grid)
!     40.03, Jul. 00: LINK1 and LINK2 in argumentlist replaced by LINK
!     40.09, Nov. 99: IMPLICIT NONE added, Method corrected
!                     Reflection option for obstacle added
!     40.14, Dec. 00: Reflection call corrected: reduced to neighbouring
!                     linepiece of obstacle (bug fix 40.11D)
!            Jan. 01: Constant waterlevel taken into account as well (bug fix 40.11E)
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added
!     40.13, Aug. 02: subroutine restructured:
!                     loop in reflection procedure changed to avoid double
!                     reflection
!                     argument list of subr REFLECT revised
!                     argument SPCDIR added
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.31, Oct. 03: changes w.r.t. obstacles
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.66, Mar. 07: extension with d'Angremond and Van der Meer transmission
!     40.80, Mar. 08: extension to unstructured grids
!
!  2. Purpose
!
!      take the value of transmission coefficient given
!      by the user in case obstacle TRANSMISSION
!
!      or
!
!      compute the transmision coeficient in case obstacle DAM
!      based on Goda (1967) [from Seelig (1979)]
!      or d'Angremond and Van der Meer formula's (1996)                   40.66
!
!      if reflections are switched on, calculate sourceterm in            40.09
!      subroutine REFLECT                                                 40.09
!
!  3. Method
!
!     Calculate transmission coefficient based on Goda (1967)             40.09
!     from Seelig (1979)                                                  40.09
!     Kt = 0.5*(1-sin {pi/(2*alpha)*(WATHIG/Hi +beta)})
!     where
!     Kt         transmission coefficient
!
!     alpha,beta coefficients dependent on structure of obstacle
!                and waves
!     WATHIG     = F = h-d is the freeboard of the dam, where h is the    40.09
!                crest level of the dam above the reference level and d   40.09
!                is the mean water level (relative to reference level)    40.09
!     Hi         incident (significant) wave height                       40.09
!                                                                         40.09
!     If reflection are switched on and obstacle is not exactly on line   40.09
!     of two neighbouring gridpoints, calculate reflections               40.09
!
!  4. Argument variables
!
!     AC2      input     Action density array                             40.09
!     ANYBIN   input     Set a particular bin TRUE or FALSE depending on  40.09
!                        SECTOR                                           40.09
!     CAX      input     Propagation velocity                             40.09
!     CAY      input     Propagation velocity                             40.09
!     CHS      input     Hs in all computational grid points
!     KCGRD    input     Grid address of points of computational stencil
!     LINK     input     indicates whether link in stencil                40.03
!                        crosses an obstacle                              40.03
!     OBREDF   output    Array of action density reduction coefficients
!                        (reduction at the obstacle)
!     REFLSO   inp/outp  contribution to the source term of action        40.41
!                        balance equation due to reflection               40.41
!     RDX,RDY  input     Array containing spatial derivative coefficients 40.09
!     WLEV2    input     Water level in grid points
!
      INTEGER  KGRPNT(MXC,MYC)
      INTEGER  LINK(2)                                                    40.03
      REAL     CHS(MCGRD), OBREDF(MDC,MSC,2), WLEV2(MCGRD)                30.70
      REAL     :: AC2(MDC,MSC,MCGRD)                                      40.09 40.22
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL     :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)                40.09 40.22
      REAL     :: REFLSO(MDC,MSC), RDX(10), RDY(10)                       40.41 40.09 40.22 40.08
      REAL     :: SPCSIG(MSC), SPCDIR(MDC,6)                              40.13 40.28
      LOGICAL  :: ANYBIN(MDC,MSC)                                         40.09
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ALOW     Lower limit for FVH
!     BK       crest width                                                40.66
!     BUPL     Upper limit for FVH
!     BVH      Bk/Hsin                                                    40.66
!     FD1      Coeff. for freq. dep. reflection: vertical displacement    40.28
!     FD2      Coeff. for freq. dep. reflection: shape parameter          40.28
!     FD3      Coeff. for freq. dep. reflection: directional coefficient  40.28
!     FD4      Coeff. for freq. dep. reflection: bending point of freq.   40.28
!     FVH      WATHIG/Hsin
!     HGT      elevation of top of obstacle above reference level
!     HSIN     incoming significant wave height
!     ID       counter in directional space
!     IENT     number of entries of this subroutine
!     ILINK    indicates which link is analyzed: 1 -> neighbour in x
!                                                2 -> neighbour in y
!     IS       counter in frequency space
!     ITRAS    indicates kind of obstacle: 0 -> constant transm
!                                          1 -> dam, Goda
!                                          2 -> dam, d'Angremond and      40.66
!                                                    Van der Meer         40.66
!     JP       counter for number of corner points of obstacles
!     L0P      wave length in deep water                                  40.66
!     LREFDIFF indicates whether reflected energy should be               40.18
!              scattered (1) or not (0)                                   40.18
!     LREFL    if LREFL=0, no reflection; if LREFL=1, constant            40.13
!              reflection coeff.                                          40.13
!     LRFRD    Indicates whether frequency dependent reflection is        40.28
!              active (#0.) or not (=0.)                                  40.28
!     NMPO     link number
!     NUMCOR   number of corner points of obstacle
!     OBET     user defined coefficient (beta) in formulation of
!              Goda/Seelig (1967/1979)
!     OBHKT    transmission coefficient in terms of wave height
!     OGAM     user defined coefficient (alpha) in formulation of
!              Goda/Seelig (1967/1979)
!     POWN     user defined power of redistribution function              40.28
!     REFLCOEF reflection coefficient in terms of action density
!     REFLTST  used to test Refl^2+Transm^2 <=1                           40.13
!     SLOPE    slope of obstacle                                          40.66
!     SQRTREF  dummy variable
!     TRCF     transmission coefficient in terms of action density
!              (user defined or calculated (in terms of waveheight))
!     X1, Y1   user coordinates of one end of grid link
!     X2, Y2   user coordinates of other end of grid link
!     X3, Y3   user coordinates of one end of obstacle side
!     X4, Y4   user coordinates of other end of obstacle side
!     XCGRID   Coordinates of computational grid in x-direction
!     XI0P     breaker parameter                                          40.66
!     XOBS     x-coordinate of obstacle point                             40.80
!     XONOBST  Indicates whether computational point (X1,Y1) is on        40.14
!              obstacle                                                   40.14
!     XV       x-coordinate of vertex of face                             40.80
!     YCGRID   Coordinates of computational grid in y-direction
!     YOBS     y-coordinate of obstacle point                             40.80
!     YV       y-coordinate of vertex of face                             40.80
!     WATHIG   freeboard of the dam (= HGT-waterlevel)
!
      INTEGER    ID, IENT, ILINK, ITRAS, IS, JP, ICGRD, LREFL,
     &           NUMCOR, NMPO, ISIGM
      INTEGER    LREFDIFF, LRFRD                                          40.31
      REAL       ALOW, BUPL, FVH, HGT, HSIN, OBET, OBHKT,
     &           SLOPE, BK, L0P, XI0P, BVH, EMAX, ETD, TP,                40.66
     &           FAC1, FAC2,
     &           POWN, OGAM, REFLCOEF,                                    40.18 40.09
     &           TRCF, X1, X2, X3, X4, Y1, Y2, Y3, Y4, WATHIG
      REAL       FD1, FD2, FD3, FD4                                       40.31 40.28
      REAL       SQRTREF                                                  40.09
      LOGICAL    XONOBST                                                  40.14
      LOGICAL :: REFLTST                                                  40.13
      REAL       XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      INTEGER    ICC, JJ
      REAL    :: XOBS(2), XV(2), YOBS(2), YV(2)                           40.80
      LOGICAL :: SwanCrossObstacle                                        40.80
      TYPE(OBSTDAT), POINTER :: COBST                                     40.31
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     REFLECT          Computes effect of reflection
!     TCROSS           Searches for crossing point if exist               40.14
!
      LOGICAL TCROSS                                                      40.14
!
!  9. Subroutines calling
!
!     SWOMPU                                                              30.70
!
! 10. Error messages
!
! 11. Remarks
!
!     Here the formulation of the transmission coefficients concerns the  40.09
!     ratio of action densities!                                          40.09
!
! 12. Structure
!
!     ------------------------------------------------------------------
!     For both links from grid point (X1,Y1) do                           40.13
!         calculate transmission coefficients                             40.09
!         assign values to OBREDF                                         40.13
!         If there is reflection                                          40.13
!         Then select obstacle side which crosses the grid link           40.13
!              calculate reflection source terms                          40.09
!     ------------------------------------------------------------------  40.13
!
! 13. Source text
! ======================================================================
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWTRCF')
!
      REFLTST = .TRUE.                                                    40.13
      DO 90 ILINK = 1 ,2
!       default transmission coefficient
        TRCF = 1.
        NMPO = LINK(ILINK)
        IF (NMPO .EQ. 0) THEN
           OBREDF(1:MDC,1:MSC,ILINK) = 1.                                 40.13
           GOTO 90
        END IF
!       incoming wave height
        HSIN = CHS(KCGRD(ILINK+1))                                        40.03
        IF (HSIN.LT.0.1E-4) HSIN = 0.1E-4
        COBST => FOBSTAC                                                  40.31
        DO JJ = 1, NMPO-1                                                 40.31
           IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT                      40.31
           COBST => COBST%NEXTOBST                                        40.31
        END DO                                                            40.31
        ITRAS  = COBST%TRTYPE                                             40.31
        IF (ITRAS .EQ. 0) THEN
!       constant transmission coefficient
!         User defined transmission coefficient concerns ratio of         40.09
!         waveheights, so:                                                40.09
          OBHKT = COBST%TRCOEF(1)                                         40.31
          TRCF  = OBHKT * OBHKT                                           40.09
        ELSE IF (ITRAS .EQ. 1) THEN
!       transmission coefficient according to Goda and Seelig
          HGT    =  COBST%TRCOEF(1)                                       40.31
          OGAM   =  COBST%TRCOEF(2)                                       40.31
          OBET   =  COBST%TRCOEF(3)                                       40.31
!         level of dam above the water surface (freeboard):
          WATHIG =  HGT - WLEV2(KCGRD(1)) - WLEV                          40.14 30.70
!
!         *** Here the transmission coeff. is that of Goda and Seelig ***
          FVH  = WATHIG/HSIN
          ALOW = -OBET-OGAM                                               40.09
          BUPL = OGAM-OBET
!
          IF (FVH.LT.ALOW) FVH = ALOW
          IF (FVH.GT.BUPL) FVH = BUPL
          OBHKT = 0.5*(1.0-SIN(PI*(FVH+OBET)/(2.0*OGAM)))
          IF (TESTFL) WRITE (PRTEST, 20) IXCGRD(1)+MXF-2,                 40.30
     &           IYCGRD(1)+MYF-2,                                         40.30
     &           ILINK, HGT, WATHIG, HSIN, OBHKT                          40.01
  20      FORMAT (' test SWTRCF ', 2X, 3I5, ' dam level=', F6.2,
     &            ' depth=', F6.2, ' Hs=', F6.2, ' transm=', F6.3)        40.01
          IF (TESTFL .AND. ITEST.GE.140) WRITE (PRTEST, 22)
     &           OGAM, OBET, ALOW, BUPL, FVH                              40.01
  22      FORMAT (8X, 5E12.4)                                             40.01
!
!         Formulation of Goda/Seelig concerns ratio of waveheights.       40.09
!         Here we use action density so:                                  40.09
          TRCF = OBHKT*OBHKT
        ELSE IF (ITRAS.EQ.2) THEN                                         40.66
!       d'Angremond and Van der Meer formulae (1996)                      40.66
          HGT   = COBST%TRCOEF(1)
          SLOPE = COBST%TRCOEF(2)
          BK    = COBST%TRCOEF(3)
!
!         level of dam above the water surface:
          WATHIG =  HGT - WLEV2(KCGRD(1)) - WLEV
!
!         compute peak frequency of incoming wave
          EMAX = 0.
          ISIGM = -1
          DO IS = 1, MSC
             ETD = 0.
             DO ID = 1, MDC
                ETD = ETD + SPCSIG(IS)*AC2(ID,IS,KCGRD(ILINK+1))*DDIR
             END DO
             IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
             END IF
          END DO
          IF (ISIGM.LE.0) ISIGM=MSC
          TP=2.*PI/SPCSIG(ISIGM)
!
!         compute breaker parameter
          L0P  = MAX(1.E-8,1.5613*TP*TP)
          XI0P = TAN(SLOPE*PI/180.)/SQRT(HSIN/L0P)
!
!         compute transmission coefficient
          FVH = WATHIG/HSIN
          BVH = BK/HSIN
          IF (BVH.EQ.0.) THEN
             OBHKT= -0.40*FVH
             IF (OBHKT.LT.0.075) OBHKT = 0.075
             IF (OBHKT.GT.0.900) OBHKT = 0.9
          ELSE IF (BVH.LT.8.) THEN
             OBHKT= -0.40*FVH + 0.64*(BVH**(-0.31))*(1.-EXP(-0.50*XI0P))
             IF (OBHKT.LT.0.075) OBHKT = 0.075
             IF (OBHKT.GT.0.900) OBHKT = 0.9
          ELSE IF (BVH.GT.12.) THEN
             OBHKT= -0.35*FVH + 0.51*(BVH**(-0.65))*(1.-EXP(-0.41*XI0P))
             IF (OBHKT.GT.0.93-0.006*BVH) OBHKT = 0.93-0.006*BVH
             IF (OBHKT.LT.0.05          ) OBHKT = 0.05
          ELSE
!            linear interpolation
             FAC1 = -0.40*FVH + 0.64*( 8.**(-0.31))*(1.-EXP(-0.50*XI0P))
             IF (FAC1.LT.0.075) FAC1 = 0.075
             IF (FAC1.GT.0.900) FAC1 = 0.9
             FAC2 = -0.35*FVH + 0.51*(12.**(-0.65))*(1.-EXP(-0.41*XI0P))
             IF (FAC2.LT.0.050) FAC2 = 0.050
             IF (FAC2.GT.0.858) FAC2 = 0.858
             OBHKT = 3.*FAC1 - 2.*FAC2 + BVH*(FAC2 - FAC1)/4.
          END IF
          IF(OPTG.NE.5) THEN
            X1 = XCGRID(IXCGRD(1),IYCGRD(1))+XOFFS
            Y1 = YCGRID(IXCGRD(1),IYCGRD(1))+YOFFS
          ELSE
            X1 = xcugrd(vs(1))+XOFFS
            Y1 = ycugrd(vs(1))+YOFFS
          ENDIF
          WRITE (PRINTF, 23) X1, Y1,
     &       ILINK, HGT, WATHIG, HSIN, SQRT(L0P/1.5613), XI0P, OBHKT      40.66
  23      FORMAT (' Transmission: ', 2X, 2F12.4, I5, ' dam level=',F6.2,
     &            ' board=', F6.2, ' Hs=', F6.2, ' Tp=', F6.2,
     &            ' Xi0p=', F6.3, ' Kt=', F6.3)                           40.66
          IF (TESTFL .AND. ITEST.GE.140) WRITE (PRTEST, 22)
     &           TAN(SLOPE*PI/180.), FVH, BVH, L0P, XI0P
!
!         Formulation of d'Angremond concerns ratio of waveheights
!         Here we use action density so:
          TRCF = OBHKT*OBHKT
        ENDIF
!       assign values to array OBREDF
        DO IS = 1, MSC
          DO ID = 1, MDC
            OBREDF(ID,IS,ILINK) = TRCF
          ENDDO
        ENDDO
!
!     *** REFLECTION ****
!     *** X1 X2 X3 ETC. are the coordinates of point according ***
!     *** with the scheme in the function TCROSS header        ***        40.04
        LREFL = COBST%RFTYP1                                              40.13
        IF ( LREFL.GT.0. ) THEN                                           40.31
!         Reflections are activated                                       40.09
          SQRTREF  = COBST%RFCOEF(1)                                      40.31
          REFLCOEF = SQRTREF * SQRTREF                                    40.09
          LREFDIFF = COBST%RFTYP2                                         40.31
          POWN     = COBST%RFCOEF(2)                                      40.31
          FD1      = COBST%RFCOEF(3)                                      40.31
          FD2      = COBST%RFCOEF(4)                                      40.31
          FD3      = COBST%RFCOEF(5)                                      40.31
          FD4      = COBST%RFCOEF(6)                                      40.31
          LRFRD    = COBST%RFTYP3                                         40.31
!
          IF (OPTG.NE.5) THEN                                             40.80
!            determine (X1,Y1) and (X2,Y2)
             ICC   = KCGRD(1)                                             40.09
             ICGRD = 0                                                    40.09
             IF ( ICC.GT.1 ) THEN                                         40.09
               X1 = XCGRID(IXCGRD(1),IYCGRD(1))                           40.09
               Y1 = YCGRID(IXCGRD(1),IYCGRD(1))                           40.09
               IF (KGRPNT(IXCGRD(ILINK+1),IYCGRD(ILINK+1)).GT.1) THEN     40.09
                 X2    = XCGRID(IXCGRD(ILINK+1),IYCGRD(ILINK+1))          40.09
                 Y2    = YCGRID(IXCGRD(ILINK+1),IYCGRD(ILINK+1))          40.09
                 ICGRD = KCGRD(ILINK+1)                                   40.09
               ENDIF
             ENDIF
             IF (ICGRD.EQ.0) GOTO 90                                      40.13
!            select obstacle side crossing the grid link
             X3 = COBST%XCRP(1)                                           40.31
             Y3 = COBST%YCRP(1)                                           40.31
             NUMCOR = COBST%NCRPTS                                        40.31
             DO JP = 2, NUMCOR                                            40.09
               X4 = COBST%XCRP(JP)                                        40.31
               Y4 = COBST%YCRP(JP)                                        40.31
               IF (TCROSS(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XONOBST)) GOTO 70       40.14
               X3 = X4                                                    40.09
               Y3 = Y4                                                    40.09
             ENDDO
          ELSE                                                            40.80
!           determine begin and end points of link                        40.80
            X1 = xcugrd(vs(1))                                            40.80
            Y1 = ycugrd(vs(1))                                            40.80
            X2 = xcugrd(vs(ILINK+1))                                      40.80
            Y2 = ycugrd(vs(ILINK+1))                                      40.80
            XV(1) = X1                                                    40.80
            YV(1) = Y1                                                    40.80
            XV(2) = X2                                                    40.80
            YV(2) = Y2                                                    40.80
!           select obstacle side crossing the grid link                   40.80
            X3 = COBST%XCRP(1)                                            40.80
            Y3 = COBST%YCRP(1)                                            40.80
            XOBS(1) = X3                                                  40.80
            YOBS(1) = Y3                                                  40.80
            DO JP = 2, COBST%NCRPTS                                       40.80
               X4 = COBST%XCRP(JP)                                        40.80
               Y4 = COBST%YCRP(JP)                                        40.80
               XOBS(2) = X4                                               40.80
               YOBS(2) = Y4                                               40.80
               IF ( SwanCrossObstacle( XV, YV, XOBS, YOBS ) ) GOTO 70     40.80
               X3 = X4                                                    40.80
               Y3 = Y4                                                    40.80
               XOBS(1) = X3                                               40.80
               YOBS(1) = Y3                                               40.80
            ENDDO                                                         40.80
          ENDIF                                                           40.80
!         no crossing found, skip procedure
          GOTO 90

  70      CALL REFLECT(AC2, REFLSO, X1, Y1, X2, Y2,                       40.13
     &                 X3, Y3, X4, Y4, CAX,                               40.13
     &                 CAY, RDX, RDY, ILINK,                              40.13
     &                 REFLCOEF, LREFDIFF, POWN, ANYBIN,                  40.13
     &                 LRFRD, SPCSIG, SPCDIR, FD1, FD2, FD3, FD4,         40.13
     &                 OBREDF, REFLTST)                                   40.13
        END IF
!
        IF (ITEST .GE. 120)  WRITE (PRTEST,10)
     &  IXCGRD(1)-1, IYCGRD(1)-1, NMPO, TRCF
  10    FORMAT(' SWTRCF: Point=', 2I5, ' NMPO  = ', I5, ' transm ',       40.03
     &  F8.3)                                                             40.03
  90  CONTINUE
      IF (.NOT.REFLTST) THEN                                              40.13
        CALL MSGERR(3,'Kt^2 + Kr^2 > 1 ')                                 40.13
        IF (ITEST.LT.50) THEN                                             40.13
          WRITE (PRTEST, 74) IXCGRD(1)-1, IYCGRD(1)-1                     40.13
  74      FORMAT (' Kt^2 + Kr^2 > 1 in grid point:', 2I4)                 40.13
        ENDIF                                                             40.13
      ENDIF                                                               40.13
      RETURN
!     * end of SUBROUTINE SWTRCF
      END
!
!************************************************************************
!                                                                       *
      SUBROUTINE REFLECT (AC2, REFLSO, X1, Y1, X2, Y2, X3, Y3,            40.13
     &                    X4, Y4, CAX, CAY, RDX, RDY,                     40.13
     &                    ILINK, REF0, LREFDIFF, POWN, ANYBIN,            40.13
     &                    LRFRD, SPCSIG, SPCDIR, FD1, FD2, FD3, FD4,      40.13
     &                    OBREDF, REFLTST)                                40.13
!                                                                       *
!************************************************************************
!                                                                         40.09
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
!                                                                         40.09
      IMPLICIT NONE                                                       40.09
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
!  0. Authors                                                             40.09
!                                                                         40.09
!     40.09  Annette Kieftenburg                                          40.09
!     40.13  Nico Booij
!     40.18  Annette Kieftenburg                                          40.18
!     40.28  Annette Kieftenburg                                          40.28
!     40.38  Annette Kieftenburg                                          40.38
!     40.41: Marcel Zijlema
!                                                                         40.09
!  1. Updates                                                             40.09
!                                                                         40.09
!     40.09, Nov. 99: Subroutine created                                  40.09
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added                40.28
!     40.38, Feb. 02: Diffuse reflection against obstacles added          40.38
!     40.13, Sep. 02: Subroutine restructured
!                     assumptions changed: reflected energy can come
!                     from Th_norm-PI/2 to Th_norm+PI/2
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent   40.08
!                     with other subroutines                              40.08
!     40.13, Nov. 03: test on refl + transm added
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!                                                                         40.09
!  2. Purpose                                                             40.09
!                                                                         40.09
!     Computation of REFLECTIONS near obstacles                           40.09
!                                                                         40.09
!  3. Method                                                              40.09
!                                                                         40.09
!     Determine the angle of the obstacle,                                40.09
!     Determine the angles between which reflections should be taken      40.09
!     into account                                                        40.09
!     Determine redistribution function                                   40.18
!     determine expression of reflection coefficient for frequency
!     dependency, if appropriate                                          40.28
!     Determine reflected action density (corrected for angle obstacle    40.09
!     and if option is on: redistribute energy)                           40.18
!     Add reflected spectrum to contribution for the right hand side      40.41 40.09
!     of matrix equation                                                  40.09
!                                                                         40.09
!  4. Modules used                                                        40.18
!                                                                         40.18
!     --                                                                  40.18
!                                                                         40.18
!  5. Argument variables                                                  40.09
!                                                                         40.09
!     AC2      inp  action density
!     ANYBIN   inp  Determines whether a bin fall within a sweep
!     CAX      inp  Propagation velocity in x-direction                   40.09
!     CAY      inp  Propagation velocity in y-direction                   40.09
!     FD1      inp  Coeff. freq. dep. reflection: vertical displacement   40.28
!     FD2      inp  Coeff. freq. dep. reflection: shape parameter         40.28
!     FD3      inp  Coeff. freq. dep. reflection: directional coefficient 40.28
!     FD4      inp  Coeff. freq. dep. reflection: bending point of freq.  40.28
!     ILINK    inp  Indicates which link is analyzed: 1 -> neighbour in x 40.09
!                                                     2 -> neighbour in y 40.09
!     LREFDIFF inp  Indicates whether reflected energy should be          40.18
!                   scattered (1) or not (0)                              40.18
!     LRFRD    inp  Indicates whether frequency dependent reflection is   40.28
!                   active (#0.) or not (=0.)                             40.28
!     OBREDF   inp  transmission coefficients
!     POWN     inp  User defined power of redistribution function         40.18
!     REF0     inp  reflection coefficient in terms of action density
!     REFLSO   i/o  contribution to the source term due to reflection     40.41
!     REFLTST  i/o  used to test Refl^2+Transm^2 <=1                      40.13
!     RDX,RDY  inp  Array containing spatial derivative coefficients      40.09
!     SPCDIR(*,1)   spectral directions (radians)
!     SPCSIG   inp  Relative frequency (= 2*PI*Freq.)                     40.28
!     X1, Y1   inp  Coordinates of computational grid point under         40.09
!                   consideration                                         40.09
!     X2, Y2   inp  Coordinates of computational grid point neighbour     40.09
!     X3, Y3   inp  User coordinates of one end of obstacle side          40.09
!     X4, Y4   inp  User coordinates of other end of obstacle side        40.09
!                                                                         40.09
      REAL       :: AC2(MDC,MSC,MCGRD)                                    40.18
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL       :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)              40.18 40.22
      REAL       :: REFLSO(MDC,MSC), OBREDF(MDC,MSC,2)                    40.41 40.18
      REAL       :: RDX(10), RDY(10)                                      40.18 40.08
      REAL       :: FD1, FD2, FD3, FD4, SPCSIG(MSC), SPCDIR(MDC,6)        40.41 40.28
      REAL       :: REF0                                                  40.18
      REAL       :: X1, X2, X3, X4, Y1, Y2, Y3, Y4                        40.18
      LOGICAL    :: ANYBIN(MDC,MSC)                                       40.18
      INTEGER    :: ILINK                                                 40.18
      REAL       :: POWN                                                  40.41 40.18
      INTEGER    :: LREFDIFF, LRFRD                                       40.41
      LOGICAL    :: REFLTST                                               40.13
!                                                                         40.18
      INTENT (IN)     AC2, CAX, CAY, OBREDF,                              40.18
     &                FD1, FD2, FD3, FD4, LRFRD,                          40.28
     &                RDX, RDY, SPCSIG, SPCDIR, X1, X2,                   40.18
     &                X3, X4, Y1, Y2, Y3, Y4, ANYBIN, ILINK,              40.18
     &                POWN, LREFDIFF                                      40.18
      INTENT (IN OUT) REF0, REFLSO                                        40.41 40.18
!                                                                         40.09
!  6. Parameter variables                                                 40.09
!                                                                         40.18
!  7. Local variables                                                     40.09
!                                                                         40.09
      REAL :: AC2REF     ! reflected action density of one spectral bin
      REAL :: BETA         ! local angle of obstacle                      40.09
      REAL :: EPS                                                         40.09
      REAL, ALLOCATABLE :: PRDIF(:)   ! scattering filter                 40.13
      REAL    :: TH_INC         ! direction of incident wave
      REAL    :: TH_NORM        ! direction of normal to obstacle
      REAL    :: TH_OUT         ! direction of outgoing wave
      REAL    :: IANG           ! angle divided by DTheta (DDIR)
      REAL    :: SUMRD          ! sum of PRDIF array
      REAL    :: W1, W2         ! interpolation coefficients

      INTEGER :: ID             ! counter of directions
      INTEGER :: IS             ! counter of frequencies
      INTEGER :: MAXIDR         ! width of scattering filter
      INTEGER :: IDR            ! relative directional counter
      INTEGER :: ID_I1, ID_I2   ! counters of incoming directions
      INTEGER :: IDA, IDB       ! counters of incoming directions
      INTEGER :: IENT           ! number of entries
!                                                                         40.18
!  8. Subroutines used                                                    40.09
!                                                                         40.09
!  9. Subroutines calling                                                 40.09
!                                                                         40.09
!     SWTRCF                                                              40.09
!                                                                         40.09
! 10. Error messages                                                      40.09
!                                                                         40.09
!     if obstacle linepiece is of length < EPS                            40.09
!                                                                         40.09
! 11. Remarks                                                             40.09
!                                                                         40.09
!    -In case the obstacle cuts exactly through computational grid point, 40.09
!     the obstacle should be moved a bit with subroutine OBSTMOVE.        40.09
!    -The length of the obstacle linepiece is assumed to be               40.09
!     'long enough' compared to grid resolution (> 0.5*sqrt(dx^2+dy^2))   40.09
!     (if this restriction is violated, the reflections due to an obsta-  40.09
!     cle of one straight line can be very different from a similar line  40.09
!     consisting of several pieces (because only the directions of the    40.09
!     spectrum that are directed towards the obstacle linepiece are       40.09
!     reflected).                                                         40.09
!    -There should be only one intersection per computational gridcell.   40.09
!     Therefore it is better to avoid sharp edges in obstacles.           40.18
!                                                                         40.09
! 12. Structure                                                           40.09
!                                                                         40.09
!     -----------------------------------------------------------------
!     Determine angle of obstacle, Beta                                   40.09
!     Determine angle of normal from (X1,Y1) to obstacle                  40.13
!     If there is constant diffuse reflection
!     Then determine scattering distribution
!     -----------------------------------------------------------------
!     For all frequencies do
!         If amount of reflection varies with frequency
!         Then determine reflection coefficient
!         -------------------------------------------------------------
!         If there is diffuse reflection
!         Then if scattering varies with frequency
!              Then determine power of cos
!              --------------------------------------------------------
!              Determine distribution
!         -------------------------------------------------------------
!         For all active directions do
!             Determine specular incoming direction
!             For directions of scattering filter do
!                 multiply incoming action with scattering coefficient
!                 add this to array AC2REF
!     -----------------------------------------------------------------
!                                                                         40.09
!     Add reflected spectrum to right hand side of matrix equation        40.09
!                                                                         40.09
! 13. Source text                                                         40.09
!                                                                         40.18
      SAVE     IENT                                                       40.09
      DATA     IENT /0/                                                   40.09
      CALL STRACE (IENT, 'REFLECT')                                       40.09
!                                                                         40.09
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))             40.09
      IF (EPS ==0) EPS = TINY(X1)                                         40.09

      IF (LREFDIFF .EQ. 0) THEN
        ALLOCATE (PRDIF(0:0))
        MAXIDR = 0
        PRDIF(0) = 1.
      ELSE
        MAXIDR = MDC/2
        ALLOCATE (PRDIF(0:MDC/2))
      ENDIF
!                                                                         40.09
!     Determine angle of obstacle BETA, and related                       40.09
!                                                                         40.09
      IF (.NOT. ((ABS(Y4-Y3).LE.EPS) .AND. (ABS(X4-X3).LE.EPS)) ) THEN    40.09
        BETA = ATAN2((Y4-Y3),(X4-X3))                                     40.09
      ELSE                                                                40.09
        CALL MSGERR (3, 'Obstacle contains line piece of length 0!')      40.09
      END IF                                                              40.09
!     determine direction of normal                         (4)
!     this is the normal from (X1,Y1)                        |
!     towards the obstacle                                   |
!     see sketch to establish sign                    (2)----+------(1)
!                                                            |
!                                                           (3)
      IF ((X1-X3)*(Y4-Y3)-(Y1-Y3)*(X4-X3) .GT. 0.) THEN                   40.13
        TH_NORM = BETA + 0.5*PI
      ELSE
        TH_NORM = BETA - 0.5*PI
      ENDIF                                                               40.13

!     prepare directional filter in case of diffuse reflection
      IF (LREFDIFF.EQ.1) THEN
        PRDIF(1:MAXIDR) = 0.
        PRDIF(0) = 1.
        SUMRD = 1.
        DO ID = 1, MAXIDR
          PRDIF(ID) = (COS(ID*DDIR))**POWN
          IF (PRDIF(ID) .GT. 0.01) THEN
            SUMRD = SUMRD + 2.*PRDIF(ID)
          ELSE
            MAXIDR = ID-1
            GOTO 23
          ENDIF
        ENDDO
  23    DO ID = 0, MAXIDR
          PRDIF(ID) = PRDIF(ID) / SUMRD
        ENDDO
        IF (TESTFL .AND. ITEST.GE.50) THEN
          WRITE (PRTEST, 27) POWN, MAXIDR
  27      FORMAT (' power scattering filter:', F4.1, I3)
          IF (ITEST.GE.130) WRITE (PRTEST, 28)
     &    (PRDIF(IDR), IDR=0, MAXIDR)
  28      FORMAT (10 F7.3)
        ENDIF
      ENDIF

      DO IS = 1, MSC
        IF (LRFRD.EQ.1) THEN
!         amount of reflection varies with wave frequency
          REF0 = FD1 +                                                    40.13
     &           FD2/PI * ATAN2(PI*FD3*(SPCSIG(IS)-FD4),FD2)              40.13
          IF (REF0 > 1.) REF0 = 1.                                        40.28
          IF (REF0 < 0.) REF0 = 0.                                        40.28
!         >>> should REF0 not be squared? <<<
        ENDIF
!       check whether reflection + transmission <= 1                      40.13
        DO ID = 1, MDC                                                    40.13
          IF ((REF0 + OBREDF(ID,IS,ILINK)) .GT. 1.) THEN                  40.13
            REFLTST = .FALSE.                                             40.13
            IF (ITEST.GE.50) THEN                                         40.13
              WRITE (PRTEST, 72) IXCGRD(1)-1, IYCGRD(1)-1, ILINK,         40.13
     &               IS, ID, REF0, OBREDF(ID,IS,ILINK)                    40.13
  72          FORMAT (' Refl+Transm>1 in ', 2I4, 2X, 3I3, 2X, 2F6.2)      40.13
            ENDIF                                                         40.13
          ENDIF                                                           40.13
        ENDDO                                                             40.13

        IF (LREFDIFF.EQ.2) THEN
!         spreading varies with frequency; not yet implemented
        ENDIF
        DO ID = 1, MDC
          IF (ANYBIN(ID,IS)) THEN
            AC2REF = 0.
            TH_OUT = SPCDIR(ID,1)
!           corresponding incident direction (assuming specular reflection)
            TH_INC = 2.*BETA-TH_OUT
!           determine counter for which direction is TH_INC:
            IANG = MOD (TH_INC-SPCDIR(1,1), 2.*PI) / DDIR
!           incident angle is between ID_I1 and ID_I2
            ID_I1 = 1 + INT (IANG)
            ID_I2 = ID_I1+1
!           W1 and W2 are weighting coefficients for the above directions
!           by linear interpolation
            W2 = IANG + 1. - REAL(ID_I1)
            W1 = 1. - W2
            DO IDR = -MAXIDR, MAXIDR
              IDA = ID_I1 + IDR
              IDB = ID_I2 + IDR
              IF (FULCIR) THEN
                IDA = 1+MOD(2*MDC+IDA-1,MDC)
!               only outgoing reflected waves, i.e. not towards obstacle
                IF (COS(TH_NORM-SPCDIR(IDA,1)) .GT. 0.)
     &          AC2REF = AC2REF +
     &          REF0 * W1 * PRDIF(ABS(IDR)) * AC2(IDA,IS,KCGRD(1))
                IDB = 1+MOD(2*MDC+IDB-1,MDC)
                IF (COS(TH_NORM-SPCDIR(IDB,1)) .GT. 0.)
     &          AC2REF = AC2REF +
     &          REF0 * W2 * PRDIF(ABS(IDR)) * AC2(IDB,IS,KCGRD(1))
              ELSE
                IF (IDA.GE.1 .AND. IDA.LE.MDC) THEN
                  IF (COS(TH_NORM-SPCDIR(IDA,1)) .GT. 0.)
     &            AC2REF = AC2REF +
     &            REF0 * W1 * PRDIF(ABS(IDR)) * AC2(IDA,IS,KCGRD(1))
                ENDIF
                IF (IDB.GE.1 .AND. IDB.LE.MDC) THEN
                  IF (COS(TH_NORM-SPCDIR(IDB,1)) .GT. 0.)
     &            AC2REF = AC2REF +
     &            REF0 * W2 * PRDIF(ABS(IDR)) * AC2(IDB,IS,KCGRD(1))
                ENDIF
              ENDIF
            ENDDO
!           add reflected energy to right hand side of matrix
            REFLSO(ID,IS) = REFLSO(ID,IS) + AC2REF *
     &      (RDX(ILINK)*CAX(ID,IS,1) + RDY(ILINK)*CAY(ID,IS,1))           40.13
          END IF                                                          40.09
        END DO                                                            40.09
      END DO                                                              40.09
!                                                                         40.18
      IF (LREFDIFF .GT. 0) DEALLOCATE (PRDIF)
      RETURN                                                              40.09
!     End of subroutine REFLECT                                           40.09
      END                                                                 40.09
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SSHAPE (ACLOC, SPCSIG, SPCDIR, FSHAPL, DSHAPL)
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
!            Roeland Ris
!            Roberto Padilla
!     30.73: Nico Booij
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.02: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!            Dec. 92: new for SWAN
!            Dec. 96: option MEAN freq. introduced see LOGPM
!     30.73, Nov. 97: revised in view of new boundary treatment
!     30.82, Sep. 98: Added error message in case of non-convergence
!     30.80, Oct. 98: correction suggested by Mauro Sclavo, and renames
!                     computation of tail added to improve accuracy
!     30.82, Oct. 98: Updated description of several variables
!     40.02, Oct. 00: Modified test write statement to avoid division by MS=0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculating of energy density at boundary point (x,y,sigma,theta)
!
!  3. Method (updated...)
!
!     see: M. Yamaguchi: Approximate expressions for integral properties
!          of the JONSWAP spectrum; Proc. JSCE, No. 345/II-1, pp. 149-152,
!          1984.
!
!     computation of mean period: see Swan system documentation
!
!  4. Argument variables
!
!   o ACLOC : Energy density at a point in space
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL    ACLOC(MDC,MSC)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
!
! i   DSHAPL: Directional distribution
! i   FSHAPL: Shape of spectrum:
!             =1; Pierson-Moskowitz spectrum
!             =2; Jonswap spectrum
!             =3; bin
!             =4; Gauss curve
!             (if >0: period is interpreted as peak per.
!              if <0: period is interpreted as mean per.)
!
      INTEGER FSHAPL, DSHAPL                                              40.00
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!     LSHAPE   absolute value of FSHAPL
!
      INTEGER  ID, IS, LSHAPE
!
!     PKPER    peak period                                                30.80
!     APSHAP   aux. var. used in computation of spectrum
!     AUX1     auxiliary variable
!     AUX2     auxiliary variable
!     AUX3     auxiliary variable
!     COEFF    coefficient for behaviour around the peak (Jonswap)
!     CPSHAP   aux. var. used in computation of spectrum
!     CTOT     total energy
!     CTOTT    total energy (used for comparison)
!     DIFPER   auxiliary variable used to select bin closest
!              to given frequency
!     MPER
!     MS       power in directional distribution
!     RA       action density
!     SALPHA
!     SF       frequency (Hz)
!     SF4      SF**4
!     SF5      SF**5
!     FPK      frequency corresponding to peak period (1/PKPER)           30.80
!     FPK4     FPK**4
!     SYF      peakedness parameter
!
      REAL     APSHAP, AUX1, AUX2, AUX3
      REAL     COEFF ,SYF   ,MPER  ,CTOT  ,CTOTT,PKPER  ,DIFPER
      REAL     MS
      REAL     RA    ,SALPHA,SF   ,SF4   ,SF5   ,FPK   ,FPK4, FAC
!
!     LOGPM    indicates whether peak or mean frequency is used
!     DVERIF   logical used in verification of incident direction
!
      LOGICAL  LOGPM, DVERIF                                              40.00
!
!     PSHAPE   coefficients of spectral distribution (see remarks)
!     SPPARM   array containing integral wave parameters (see remarks)
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
! 10. Error messages
!
! 11. Remarks
!
!     PSHAPE(1): SY0, peak enhancement factor (gamma) in Jonswap spectrum
!     PSHAPE(2): spectral width in case of Gauss spectrum in rad/s
!
!     SPPARM    real     input    incident wave parameters (Hs, Period,
!                                 direction, Ms (dir. spread))
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by the user (either peak or mean)      30.80
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!
!     ---------------------------------------------------------------------
!
!     In the case of a JONSWAP spectrum the initial conditions are given by
!                   _               _       _       _       _
!                  |       _   _ -4  |     |       | S - S   |
!             2    |      |  S  |    |     |       |      p  |
!          a g     |      |  _  |    |  exp|-1/2 * |________ |* 2/pi COS(T-T  )
! E(S,D )= ___  exp|-5/4 *|  S  |    | G   |       | e * S   |              wi
!      wa    5     |      |   p |    |     |_      |_     p _|
!           S      |      |_   _|    |
!                  |_               _|
!
!   where
!         S   : rel. frequency
!
!         D   : Dir. of wave component
!          wa
!
!         a   : equili. range const. (Phillips' constant)
!         g   : gravity acceleration
!
!         S   : Peak frequency
!          p
!
!         G   : Peak enhancement factor
!         e   : Peak width
!
!         T   : local wind direction
!          wi
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       case shape
!       =1:   calculate value of Pierson-Moskowitz spectrum
!       =2:   calculate value of Jonswap spectrum
!       =3:   calculate value of bin spectrum
!       =4:   calculate value of Gauss spectrum
!       else: Give error message because of wrong shape
!       ----------------------------------------------------------------
!       if LOGPM is True
!       then calculate average period
!            if it differs from given average period
!            then recalculate peak period
!                 restart procedure to compute spectral shape
!       ----------------------------------------------------------------
!       for all spectral bins do
!            multiply all action densities by directional distribution
!       ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE(IENT,'SSHAPE')
!
      IF (ITEST.GE.80) WRITE (PRTEST, 8) FSHAPL, DSHAPL,
     &      (SPPARM(JJ), JJ = 1,4)
   8  FORMAT (' entry SSHAPE ', 2I3, 4E12.4)
      IF (FSHAPL.LT.0) THEN
        LSHAPE = - FSHAPL
        LOGPM  = .FALSE.
      ELSE
        LSHAPE = FSHAPL
        LOGPM  = .TRUE.
      ENDIF
!
      IF (SPPARM(1).LE.0.)                                                40.31
     &   CALL MSGERR(1,' sign. wave height at boundary is not positive')  40.31
!
      PKPER = SPPARM(2)
      ITPER = 0
      IF (LSHAPE.EQ.3) THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF (ABS(PKPER - PI2/SPCSIG(IS)) .LT. DIFPER) THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2/SPCSIG(IS))
          ENDIF
        ENDDO
      ENDIF
!
!     compute spectral shape using peak period PKPER                      30.80
!
      FAC  = 1.
 100  FPK  = (1./PKPER)                                                   30.80
      FPK4 = FPK**4
      IF (LSHAPE.EQ.1) THEN
        SALPHA = ((SPPARM(1) ** 2) * (FPK4)) * 5. / 16.
      ELSE IF (LSHAPE.EQ.2) THEN
!       *** SALPHA = alpha*(grav**2)/(2.*pi)**4)
        SALPHA = (SPPARM(1)**2 * FPK4) /
     &             ((0.06533*(PSHAPE(1)**0.8015)+0.13467)*16.)
      ELSE IF (LSHAPE.EQ.4) THEN
        AUX1 = SPPARM(1)**2 / ( 16.* SQRT (PI2) * PSHAPE(2))
        AUX3 = 2. * PSHAPE(2)**2
      ENDIF
!
      CTOTT = 0.
      DO 300 IS = 1, MSC                                                  30.80
!
        IF (LSHAPE.EQ.1) THEN
!         *** LSHAPE = 1 : Pierson and Moskowitz ***
          SF = SPCSIG(IS) / PI2
          SF4 = SF**4
          SF5 = SF**5
          RA = (SALPHA/SF5)*EXP(-(5.*FPK4)/(4.*SF4))/(PI2*SPCSIG(IS))
          ACLOC(MDC,IS) = RA
        ELSE IF (LSHAPE.EQ.2) THEN
!         *** LSHAPE = 2 : JONSWAP ***
          SF = SPCSIG(IS)/(PI2)
          SF4 = SF**4
          SF5 = SF**5
          CPSHAP = 1.25 * FPK4 / SF4
          IF (CPSHAP.GT.10.) THEN                                         30.50
            RA = 0.
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF
          IF (SF .LT. FPK) THEN
            COEFF = 0.07
          ELSE
            COEFF = 0.09
          ENDIF
          APSHAP =  0.5 * ((SF-FPK) / (COEFF*FPK)) **2
          IF (APSHAP.GT.10.) THEN                                         30.50
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = PSHAPE(1)**PPSHAP
          ENDIF
          RA = SYF*RA/(SPCSIG(IS)*PI2)
          ACLOC(MDC,IS) = RA
          IF (ITEST.GE.120) WRITE (PRTEST, 112)
     &                 SF, SALPHA, CPSHAP, APSHAP, SYF, RA
 112      FORMAT (' SSHAPE freq. ', 8E12.4)
        ELSE IF (LSHAPE.EQ.3) THEN
!
!         *** all energy concentrated in one BIN ***
!
          IF (IS.EQ.ISP) THEN
            ACLOC(MDC,IS) = ( SPPARM(1)**2 ) /
     &                     ( 16. * SPCSIG(IS)**2 * FRINTF )
          ELSE
            ACLOC(MDC,IS) = 0.
          ENDIF
        ELSE IF (LSHAPE.EQ.4) THEN
!
!         *** energy Gaussian distributed (wave-current tests) ***
!
          AUX2 = ( SPCSIG(IS) - ( PI2 / PKPER ) )**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPCSIG(IS)
          ACLOC(MDC,IS) = RA
        ELSE
          IF (IS.EQ.1) THEN
            CALL MSGERR (2,'Wrong type for frequency shape')
            WRITE (PRINTF, *) ' -> ', FSHAPL, LSHAPE
          ENDIF
        ENDIF
        IF (ITEST.GE.10)
     &        CTOTT = CTOTT + FRINTF * ACLOC(MDC,IS) * SPCSIG(IS)**2
 300  CONTINUE
      IF (ITEST.GE.10) THEN
        IF (SPPARM(1).GT.0.01) THEN
          HSTMP = 4. * SQRT(CTOTT)
          IF (ABS(HSTMP-SPPARM(1)) .GT. 0.1*SPPARM(1))
     &    WRITE (PRINTF, 303) SPPARM(1), HSTMP
 303      FORMAT (' SSHAPE, deviation in Hs, should be ', F8.3,
     &            ', calculated ', F8.3)
        ENDIF
      ENDIF
!
!     if mean frequency was given recalculate PKPER and restart
!
      IF (.NOT.LOGPM .AND. ITPER.LT.10) THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(MDC,IS) * (SPCSIG(IS))**2
          AS3 = AS2 * SPCSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PWTAIL(1) - 1.                                           30.80
        APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              30.80
        AM0 = AM0 * FRINTF + APTAIL * AS2                                 30.80
        PPTAIL = PWTAIL(1) - 2.                                           30.80
        EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              30.80
        AM1 = AM1 * FRINTF + EPTAIL * AS3                                 30.80
!       Mean period:
        IF ( AM1.NE.0. ) THEN                                             40.31
           MPER = PI2 * AM0 / AM1
        ELSE                                                              40.31
           CALL MSGERR(3, ' first moment is zero in calculating the')     40.31
           CALL MSGERR(3, ' spectrum at boundary using param. bc.')       40.31
        END IF                                                            40.31
        IF (ITEST.GE.80) WRITE (PRTEST, 72) ITPER, SPPARM(2), MPER,
     &          PKPER
  72    FORMAT (' SSHAPE iter=', I2, '  period values:', 3F7.2)
        IF (ABS(MPER-SPPARM(2)) .GT. 0.01*SPPARM(2)) THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPARM(2) / MPER) * PKPER                              30.80
          GOTO 100
        ENDIF
      ENDIF
!
      IF (ITPER.GE.10) THEN
        CALL MSGERR(3, 'No convergence calculating the spectrum')         30.82
        CALL MSGERR(3, 'at the boundary using parametric bound. cond.')   30.82
      ENDIF
!
!     now introduce distribution over directions
!
      ADIR = PI * DEGCNV(SPPARM(3)) / 180.                                40.00
      IF (DSHAPL.EQ.1) THEN
        DSPR = PI * SPPARM(4) / 180.
        MS = MAX (DSPR**(-2) - 2., 1.)
      ELSE
        MS = SPPARM(4)
      ENDIF
      IF (MS.LT.12.) THEN
        CTOT = (2.**MS) * (GAMMA(0.5*MS+1.))**2 / (PI * GAMMA(MS+1.))
      ELSE
        CTOT =  SQRT (0.5*MS/PI) / (1. - 0.25/MS)
      ENDIF
      IF (ITEST.GE.100) THEN
        ESOM = 0.
        DO IS = 1, MSC
          ESOM = ESOM + FRINTF * SPCSIG(IS)**2 * ACLOC(MDC,IS)
        ENDDO
        GAM1 = GAMMA(0.5*MS+1.)                                                40.84
        GAM2 = GAMMA(MS+1.)                                                    40.84
        WRITE (PRTEST, *) ' SSHAPE dir ', 4.*SQRT(ABS(ESOM)),
     &        SPPARM(1), CTOT, MS, GAM1, GAM2, CTOT                            40.84 40.02
      ENDIF
      DVERIF = .FALSE.
      CTOTT = 0.
      DO ID = 1, MDC
        ACOS = COS(SPCDIR(ID,1) - ADIR)
        IF (ACOS .GT. 0.) THEN
          CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
          IF (.NOT.FULCIR) THEN
            IF (ACOS .GE. COS(DDIR)) DVERIF = .TRUE.
          ENDIF
        ELSE
          CDIR = 0.
        ENDIF
        IF (ITEST.GE.10) CTOTT = CTOTT + CDIR * DDIR
        IF (ITEST.GE.100) WRITE (PRTEST, 360) ID,SPCDIR(ID,1),CDIR
 360    FORMAT (' ID Spcdir Cdir: ',I3,3(1X,E10.4))
        DO IS = 1, MSC
          ACLOC(ID,IS) = CDIR * ACLOC(MDC,IS)
        ENDDO
      ENDDO
      IF (ITEST.GE.10) THEN
        IF (ABS(CTOTT-1.) .GT. 0.1) WRITE (PRINTF, 363) CTOTT
 363    FORMAT (' SSHAPE, integral of Cdir is not 1, but:', F6.3)
      ENDIF
      IF (.NOT.FULCIR .AND. .NOT.DVERIF)
     &   CALL MSGERR (1, 'incident direction is outside sector')
!
      RETURN
!
! End of subroutine SSHAPE
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE SINTRP (W1, W2, FL1, FL2, FL, SPCDIR, SPCSIG)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
!
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  Weimin Luo, Roeland Ris, Nico Booij            |
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
!     30.82: IJsbrand Haagsma
!     40.00: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.01, Jan. 96: New subroutine for SWAN Ver. 30.01
!     30.73, Nov. 97: revised
!     40.00, Apr. 98: procedure to maintain peakedness introduced
!     30.82, Oct. 98: Update description of several variables
!     30.82, Oct. 98: Made arguments in ATAN2 double precision to prevent
!                     underflows
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     interpolation of spectra
!
!  3. Method (updated...)
!
!     linear interpolation with peakedness maintained
!     interpolated average direction and frequency are determined
!     average direction and frequency of interpolated spectrum are determ.
!     shifts in frequency and direction are determined from spectrum 1 and
!     2 to the interpolated spectrum
!     bilinear interpolation in spectral space is used to calculate
!     contributions from spectrum 1 and 2.
!     in full circle cases interpolation crosses the boundary 0-360 degr.
!
!  4. Argument variables
!
!   o FL    : Interpolated spectrum.
! i   FL1   : Input spectrum 1.
! i   FL2   : Input spectrum 2.
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
! i   W1    : Weighting coefficient for spectrum 1.
! i   W2    : Weighting coefficient for spectrum 2.
!
      REAL    FL1(MDC,MSC), FL2(MDC,MSC), FL(MDC, MSC)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.82
      REAL    W1, W2
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!
      INTEGER  ID, IS
!
!     DOADD    indicates whether or not values have to be added
!
      LOGICAL  DOADD
!
!     ATOT1    integral over spectrum 1
!     ATOT2    integral over spectrum 2
!     AXTOT1   integral over x-component of spectrum 1
!     AXTOT2   integral over x-component of spectrum 2
!     AYTOT1   integral over y-component of spectrum 1
!     AYTOT2   integral over y-component of spectrum 2
!     ASTOT1   integral over Sigma * spectrum 1
!     ASTOT2   integral over Sigma * spectrum 2
!     ASIG1    average Sigma of spectrum 1
!     ASIG2    average Sigma of spectrum 2
!     DELD1    difference in direction between spectrum 1 and
!              the interpolated spectrum in number of directional steps
!     DELD2    same for spectrum 2
!     DELSG1   shift in frequency between spectrum 1 and interpolated
!              spectrum in number of frequency steps
!     DELSG2   same for spectrum 2
!
      REAL     ATOT1,  ATOT2,  AXTOT1, AXTOT2, AYTOT1, AYTOT2,
     &         ASTOT1, ASTOT2
      REAL     ASIG1,  ASIG2
      REAL     DELD1,  DELD2,  DELSG1, DELSG2
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!      SNEXTI, RBFILE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!      -----------------------------------------------------------------
!      If W1 close to 1
!      Then copy FL from FL1
!      Else If W2 close to 1
!           Then copy FL from FL2
!           Else determine total energy in FL1 and FL2
!                If energy of FL1 = 0
!                Then make FL = W2 * FL2
!                Else If energy of FL2 = 0
!                     Then make FL = W1 * FL1
!                     Else determine average direction of FL1 and FL2
!                          make ADIR = W1 * ADIR1 + W2 * ADIR2
!                          determine average frequency of FL1 and FL2
!                          make ASIG = W1 * ASIG1 + W2 * ASIG2
!                          determine directional shift from FL1
!                          determine directional shift from FL2
!                          determine frequency shift from FL1
!                          determine frequency shift from FL2
!                          For all spectral components do
!                              compose FL from components of FL1 and FL2
!      -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINTRP')
!
!     interpolation of spectra
!     ------------------------
!
      IF (W1.GT.0.99) THEN
        DO 101 ID=1,MDC
          DO 102 IS=1,MSC
            FL(ID,IS) = FL1(ID,IS)
 102      CONTINUE
 101    CONTINUE
      ELSE IF (W1.LT.0.01) THEN
        DO 201 ID=1,MDC
          DO 202 IS=1,MSC
            FL(ID,IS) = FL2(ID,IS)
 202      CONTINUE
 201    CONTINUE
      ELSE
        ATOT1  = 0.
        ATOT2  = 0.
        AXTOT1 = 0.
        AXTOT2 = 0.
        AYTOT1 = 0.
        AYTOT2 = 0.
        ASTOT1 = 0.
        ASTOT2 = 0.
        DO 301 ID=1,MDC
          DO 302 IS=1,MSC
            ATOT1  = ATOT1  + FL1(ID,IS)
            AXTOT1 = AXTOT1 + FL1(ID,IS) * SPCDIR(ID,2)
            AYTOT1 = AYTOT1 + FL1(ID,IS) * SPCDIR(ID,3)
            ASTOT1 = ASTOT1 + FL1(ID,IS) * SPCSIG(IS)
            ATOT2  = ATOT2  + FL2(ID,IS)
            AXTOT2 = AXTOT2 + FL2(ID,IS) * SPCDIR(ID,2)
            AYTOT2 = AYTOT2 + FL2(ID,IS) * SPCDIR(ID,3)
            ASTOT2 = ASTOT2 + FL2(ID,IS) * SPCSIG(IS)
 302      CONTINUE
 301    CONTINUE
        IF (ATOT1.LT.1.E-9) THEN
          DO 401 ID=1,MDC
            DO 402 IS=1,MSC
              FL(ID,IS) = W2*FL2(ID,IS)
 402        CONTINUE
 401      CONTINUE
        ELSE IF (ATOT2.LT.1.E-9) THEN
          DO 501 ID=1,MDC
            DO 502 IS=1,MSC
              FL(ID,IS) = W1*FL1(ID,IS)
 502        CONTINUE
 501      CONTINUE
        ELSE
!         determine interpolation factors in Theta space
          AXTOT  = W1 * AXTOT1 + W2 * AXTOT2
          AYTOT  = W1 * AYTOT1 + W2 * AYTOT2
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 509)  ATOT1, ATOT2,                            40.02
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2
 509        FORMAT (' SINTRP factors ', 8E11.4, /, 15X, 4F7.3)
          ENDIF
!         DELD1 is the difference in direction between spectrum 1 and
!         the interpolated spectrum in number of directional steps
          DELD1  = REAL(ATAN2(DBLE(AXTOT*AYTOT1 - AYTOT*AXTOT1),          30.82
     &                        DBLE(AXTOT*AXTOT1 + AYTOT*AYTOT1))) / DDIR  30.82
!         DELD2 is the difference between spectrum 2 and
!         the interpolated spectrum
          DELD2  = REAL(ATAN2(DBLE(AXTOT*AYTOT2 - AYTOT*AXTOT2),          30.82
     &                        DBLE(AXTOT*AXTOT2 + AYTOT*AYTOT2))) / DDIR  30.82
          IDD1A  = NINT(DELD1)
          RDD1B  = DELD1 - REAL(IDD1A)
          IF (RDD1B .LT. 0.) THEN
            IDD1A = IDD1A - 1
            RDD1B = RDD1B + 1.
          ENDIF
          IDD1B  = IDD1A + 1
          RDD1B  = W1 * RDD1B
          RDD1A  = W1 - RDD1B
          IDD2A  = NINT(DELD2)
          RDD2B  = DELD2 - REAL(IDD2A)
          IF (RDD2B .LT. 0.) THEN
            IDD2A = IDD2A - 1
            RDD2B = RDD2B + 1.
          ENDIF
          IDD2B  = IDD2A + 1
          RDD2B  = W2 * RDD2B
          RDD2A  = W2 - RDD2B
!
!         determine interpolation factors in Sigma space
          ASIG1  = ASTOT1 / ATOT1
          ASIG2  = ASTOT2 / ATOT2
          ATOT   = W1 * ATOT1  + W2 * ATOT2
          ASTOT  = W1 * ASTOT1 + W2 * ASTOT2
          ASIG   = ASTOT / ATOT
!
!         DELSG1 is shift in frequency between spectrum 1 and interpolated
!         spectrum in number of frequency steps
          DELSG1 = ALOG (ASIG1 / ASIG) / FRINTF
          IDS1A  = NINT(DELSG1)
          RDS1B  = DELSG1 - REAL(IDS1A)
          IF (RDS1B .LT. 0.) THEN
            IDS1A = IDS1A - 1
            RDS1B = RDS1B + 1.
          ENDIF
          IDS1B  = IDS1A + 1
          RDS1A  = 1. - RDS1B
!
!         DELSG2 is shift in frequency between spectrum 2 and interpolated
!         spectrum in number of frequency steps
          DELSG2 = ALOG (ASIG2 / ASIG) / FRINTF
          IDS2A  = NINT(DELSG2)
          RDS2B  = DELSG2 - REAL(IDS2A)
          IF (RDS2B .LT. 0.) THEN
            IDS2A = IDS2A - 1
            RDS2B = RDS2B + 1.
          ENDIF
          IDS2B  = IDS2A + 1
          RDS2A  = 1. - RDS2B
!         test output
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 510) ATOT, ATOT1, ATOT2,
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2,
     &            DELD1, DELD2, DELSG1, DELSG2
 510        FORMAT (' SINTRP factors ', 9E11.4, /, 15X, 4F7.3)
            WRITE (PRTEST, 512) IDS1A, RDS1A, IDS1B, RDS1B,
     &            IDS2A, RDS2A, IDS2B, RDS2B,
     &            IDD1A, RDD1A, IDD1B, RDD1B,
     &            IDD2A, RDD2A, IDD2B, RDD2B
 512        FORMAT (' SINTRP ', 8(I2, F7.3))
          ENDIF
!
          DO 601 ID=1,MDC
            DO 602 IS=1,MSC
              FL(ID,IS) = 0.
 602        CONTINUE
 601      CONTINUE
          DO 611 ID=1,MDC
            DOADD = .TRUE.
            ID1A = ID + IDD1A
            IF (FULCIR) THEN
              IF (ID1A.LT.1)   ID1A = ID1A + MDC
              IF (ID1A.GT.MDC) ID1A = ID1A - MDC
            ELSE
              IF (ID1A.LT.1)   DOADD = .FALSE.
              IF (ID1A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 612 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1A * RDS1A * FL1(ID1A,IS+IDS1A)
 612          CONTINUE
              DO 613 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1A * RDS1B * FL1(ID1A,IS+IDS1B)
 613          CONTINUE
            ENDIF
 611      CONTINUE
          DO 621 ID=1,MDC
            DOADD = .TRUE.
            ID1B = ID + IDD1B
            IF (FULCIR) THEN
              IF (ID1B.LT.1)   ID1B = ID1B + MDC
              IF (ID1B.GT.MDC) ID1B = ID1B - MDC
            ELSE
              IF (ID1B.LT.1)   DOADD = .FALSE.
              IF (ID1B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 622 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1B * RDS1A * FL1(ID1B,IS+IDS1A)
 622          CONTINUE
              DO 623 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD1B * RDS1B * FL1(ID1B,IS+IDS1B)
 623          CONTINUE
            ENDIF
 621      CONTINUE
          DO 631 ID=1,MDC
            DOADD = .TRUE.
            ID2A = ID + IDD2A
            IF (FULCIR) THEN
              IF (ID2A.LT.1)   ID2A = ID2A + MDC
              IF (ID2A.GT.MDC) ID2A = ID2A - MDC
            ELSE
              IF (ID2A.LT.1)   DOADD = .FALSE.
              IF (ID2A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 632 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2A * RDS2A * FL2(ID2A,IS+IDS2A)
 632          CONTINUE
              DO 633 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2A * RDS2B * FL2(ID2A,IS+IDS2B)
 633          CONTINUE
            ENDIF
 631      CONTINUE
          DO 641 ID=1,MDC
            DOADD = .TRUE.
            ID2B = ID + IDD2B
            IF (FULCIR) THEN
              IF (ID2B.LT.1)   ID2B = ID2B + MDC
              IF (ID2B.GT.MDC) ID2B = ID2B - MDC
            ELSE
              IF (ID2B.LT.1)   DOADD = .FALSE.
              IF (ID2B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 642 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2B * RDS2A * FL2(ID2B,IS+IDS2A)
 642          CONTINUE
              DO 643 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) +
     &                      RDD2B * RDS2B * FL2(ID2B,IS+IDS2B)
 643          CONTINUE
            ENDIF
 641      CONTINUE
        ENDIF
      ENDIF
!
!     Test output
      IF (ITEST.GE.80) THEN
        A1 = 0.
        A2 = 0.
        AA = 0.
        DO 801 ID=1,MDC
          DO 802 IS=1,MSC
            A1 = MAX(A1,FL1(ID,IS))
            A2 = MAX(A2,FL2(ID,IS))
            AA = MAX(AA,FL(ID,IS))
 802      CONTINUE
 801    CONTINUE
        WRITE (PRTEST, *) ' SINTRP, maxima ', A1, A2, AA
      ENDIF
!
      RETURN
!  end of subroutine of SINTRP
      END
!***********************************************************************
!                                                                      *
      REAL FUNCTION DEGCNV (DEGREE)
!                                                                      *
!***********************************************************************
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees from nautical to cartesian or vice versa.
!
!  3. METHOD
!
!       DEGCNV = 180 + dnorth - degree
!
!  4. PARAMETERLIST
!
!       DEGCNV      direction in cartesian or nautical degrees.
!       DEGREE      direction in nautical or cartesian degrees.
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!           Nautical convention           Cartesian convention
!
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
!  9. STRUCTURE
!
!     ---------------------------------
!     IF (NAUTICAL DEGREES) THEN
!       CONVERT DEGREES
!     IF (DEGREES > 360 OR < 0) THEN
!       CORRECT DEGREES WITHIN 0 - 360
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'DEGCNV')
!
      IF ( BNAUT ) THEN
        DEGCNV = 180. + DNORTH - DEGREE
      ELSE
        DEGCNV = DEGREE
      ENDIF
!
      IF (DEGCNV .GE. 360.) THEN
        DEGCNV = MOD (DEGCNV, 360.)
      ELSE IF (DEGCNV .LT. 0.) THEN
        DEGCNV = MOD (DEGCNV, 360.) + 360.
      ELSE
!       DEGCNV between 0 and 360; do nothing
      ENDIF
!
!
!     *** end of subroutine DEGCNV ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGRAD (DEGREE)
!                                                                      *
!***********************************************************************
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees to radians
!
!  3. METHOD
!
!       ANGRAD = DEGREE * PI / 180
!
!  4. PARAMETERLIST
!
!       ANGRAD      radians
!       DEGREE      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[radian] = ANGLE[degrees} * PI / 180
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGRAD')
!
      ANGRAD = DEGREE * PI / 180.
!
!
!     *** end of subroutine ANGRAD ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGDEG (RADIAN)
!                                                                      *
!***********************************************************************
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform radians to degrees
!
!  3. METHOD
!
!       ANGDEG = RADIAN * 180 / PI
!
!  4. PARAMETERLIST
!
!       RADIAN      radians
!       ANGDEG      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[degrees] = ANGLE[radians} * 180 / PI
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGDEG')
!
      ANGDEG = RADIAN * 180. / PI
!
!
!     *** end of subroutine ANGDEG ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      SUBROUTINE HSOBND (AC2   ,SPCSIG,HSIBC ,KGRPNT)                     40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_PARALL                                                        40.31
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
!     32.01: Roeland Ris
!     30.70: Nico Booij
!     40.00: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     32.01, Sep. 97: new for SWAN
!     30.72, Jan. 98: Changed number of elements for HSI to MCGRD
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: structure scheme corrected
!     40.00, Mar. 98: integration method changed (as in SNEXTI)
!                     structure corrected
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Compare computed significant wave height with the value of
!     the significant wave height as predescribed by the user. If
!     the values differ more than e.g. 10 % give an error message
!     and the gridpoints where the error has been located
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: input  Relative frequencies in computational domain in      30.72
!                    sigma-space                                          30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!       REALS:
!       ------
!       AC2        action density
!       HSI        significant wave height at boundary (using SWAN
!                  resolution (has thus not to be equal to the WAVEC
!                  significant wave height )
!       ETOT       total energy in a gridpoint
!       DS         increment in frequency space
!       DDIR       increment in directional space
!       HSC        computed wave height after SWAN computation
!       EFTAIL     contribution of tail to spectrum
!
!       INTEGERS:
!       ---------
!       KGRPNT     values of grid indices
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       TRACE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ------------------------------------------------------------------
!     for all computational grid points do                                30.70
!         if HSI is non-zero
!         then compute Hs from action density array
!              if relative difference is large than HSRERR
!              then write error message
!    -------------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      REAL      AC2(MDC,MSC,MCGRD) ,HSIBC(MCGRD)                          30.72
!
      REAL      ETOT  ,HSC                                                40.00
!
      INTEGER   ID    ,IS     ,IX     ,IY    ,INDX
!
      LOGICAL   HSRR
!
      INTEGER   KGRPNT(MXC,MYC)
!
      SAVE IENT, HSRR
      DATA IENT/0/, HSRR/.TRUE./
      CALL STRACE (IENT, 'HSOBND')
!
!     *** initializing ***
!
      HSRR = .TRUE.                                                       40.03
!
      DO IY = MYC, 1, -1
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF ( HSIBC(INDX) .GT. 1.E-25 ) THEN
!           *** compute Hs for boundary point (without tail) ***
            ETOT  = 0.
            DO ID = 1, MDC
              DO IS = 1, MSC                                              40.00
                ETOT = ETOT + SPCSIG(IS)**2 * AC2(ID,IS,INDX)             40.00
              ENDDO
            ENDDO
            IF (ETOT .GT. 1.E-8) THEN
              HSC = 4. * SQRT(ETOT*FRINTF*DDIR)
            ELSE
              HSC = 0.
            ENDIF
            HSREL = ABS(HSIBC(INDX) - HSC) / HSIBC(INDX)                  40.51 40.00
            IF (HSREL .GT. HSRERR) THEN
              IF ( HSRR ) THEN
                WRITE (PRINTF,*) ' ** WARNING : ',
     &             'Differences in wave height at the boundary'
                WRITE (PRINTF,802) HSRERR
 802            FORMAT (' Relative difference between input and ',
     &          'computation >= ', F6.2)
                WRITE (PRINTF,*) '                        Hs[m]',
     &                           '      Hs[m]      Hs[-]'
                WRITE (PRINTF,*) '    ix    iy  index   (input)',
     &                           ' (computed) (relative)'
                WRITE (PRINTF,*) ' ----------------------------',
     &                           '----------------------'
                HSRR = .FALSE.
              ENDIF
              WRITE (PRINTF,'(2(1x,I5),I7,3(1x,F10.2))')
     &                 IX+MXF-1, IY+MYF-1, INDX, HSIBC(INDX), HSC, HSREL
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      WRITE(PRINTF,*)
!
      IF ( ITEST .GE. 150 ) THEN
        WRITE(PRINTF,*) 'Values of wave height at boundary (HSOBND)'
        WRITE(PRINTF,*) '------------------------------------------'
        DO IY = MYC, 1, -1
          WRITE (PRINTF,'(13F8.3)') ( HSIBC(KGRPNT(IX,IY)), IX=1 , MXC)
        ENDDO
      ENDIF
!
!     *** end of subroutine HSOBND ***
!
      RETURN
      END
!
!*****************************************************************
!                                                                *
      SUBROUTINE CHGBAS (X1, X2, PERIOD, Y1, Y2, N1, N2,
     &                   ITEST, PRTEST)
!                                                                *
!*****************************************************************
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  G. van Vledder, N. Booij                       |
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
!  0. Update history
!
!       ver 20.48: also accomodates periodic variables such as directions
!
!  1. Purpose
!
!       change x-basis of a discretized y-function
!
!  2. Method
!
!     A piecewise constant representation of the functions is assumed
!
!     first boundaries of a cell in X1 are determined
!     then it is determined whether there are overlaps with cells
!     in X2. if so Y1*common length is added to Y2
!     Finally Y2 values are divided by cell lengths
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
!     X1       i    ra   x-coordinates of input grid
!     X2       i    ra   x-coordinates of output grid
!     PERIOD   i    r    period, i.e. x-axis is periodic if period>0
!                        e.g. spectral directions
!     Y1       i    ra   function values of input grid
!     Y2       o    ra   function values of output grid
!     N1       i    i    number of x-values of input grid
!     N2       i    i    number of x-values of output grid
!
!  4. Subroutines used
!
!     ---
!
!  5. Error messages
!
!  6. Remarks
!
!       Cell boundaries in X1 are: X1A and X1B
!       X2 is assumed to be monotonically increasing; this is checked
!       X1 is assumed to be monotonous but not necessarily increasing
!
!  7. Structure
!
!       ------------------------------------------------------------------
!       Make all values of Y2 = 0
!       For each cell in X1 do
!           determine boundaries of cell in X1
!           --------------------------------------------------------------
!           For each cell in X2 do
!               determine overlap with cell in X1; limits: RLOW and RUPP
!               add to Y2: Y1 * length of overlapping interval
!       ------------------------------------------------------------------
!       For each cell in X2 do
!           divide Y2 value by cell length
!       ------------------------------------------------------------------
!
!  8. Source text
!
      INTEGER  I1, I2, N1, N2, ITEST, PRTEST
      REAL     X1(N1), Y1(N1), X2(N2), Y2(N2), PERIOD
      REAL     X1A, X1B, X2A, X2B, RLOW, RUPP
      LOGICAL  TWICE
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'CHGBAS')
!
!     initialize output data
!
      DO I2 = 1, N2
        Y2(I2) = 0.
      ENDDO
      DO I2 = 2, N2
        IF (X2(I2).LE.X2(I2-1))
     &    CALL MSGERR (2, 'subr. CHGBAS: values of X2 not increasing')
      ENDDO
!     boundaries of the range in X2
      X2LO  = 1.5 * X2(1)  - 0.5 * X2(2)
      X2HI  = 1.5 * X2(N2) - 0.5 * X2(N2-1)
      TWICE = .FALSE.
!
!     loop over cells in X1
!
      DO 300 I1 = 1, N1
        IF (ABS(Y1(I1)) .LT. 1.E-20) GOTO 300
!
!       determine cell boundaries in X1
!
        IF (I1.EQ.1) THEN
          X1A = 1.5 * X1(1) - 0.5 * X1(2)
        ELSE
          X1A = 0.5 * (X1(I1) + X1(I1-1))
        ENDIF

        IF (I1.EQ.N1) THEN
          X1B = 1.5 * X1(N1) - 0.5 * X1(N1-1)
        ELSE
          X1B = 0.5 * (X1(I1) + X1(I1+1))
        ENDIF
!
!       swap X1A and X1B if X1A > X1B
!
        IF (X1A.GT.X1B) THEN
          RR  = X1A
          X1A = X1B
          X1B = RR
        ENDIF

        IF (PERIOD.LE.0.) THEN
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
        ELSE
!         X is periodic; move interval in X1 if necessary
          TWICE = .FALSE.
          IADD = 0
  60      IF (X1B.GT.X2HI) THEN
            X1A = X1A - PERIOD
            X1B = X1B - PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &         CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 60
          ENDIF
  70      IF (X1A.LT.X2LO) THEN
            X1A = X1A + PERIOD
            X1B = X1B + PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &           CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 70
          ENDIF
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
          IF (X1A.LT.X2LO .AND. X1A+PERIOD.LT.X2HI) TWICE = .TRUE.
          IF (X1B.GT.X2HI .AND. X1B-PERIOD.GT.X2LO) TWICE = .TRUE.
        ENDIF
!
!       loop over cells in X2
!
 100    DO 200 I2 = 1, N2

          IF (I2.EQ.1) THEN
            X2A = X2LO
          ELSE
            X2A = 0.5 * (X2(I2) + X2(I2-1))
          ENDIF

          IF (I2.EQ.N2) THEN
            X2B = X2HI
          ELSE
            X2B = 0.5 * (X2(I2) + X2(I2+1))
          ENDIF
!
!         (RLOW,RUPP) is overlapping interval of (X1A,X1B) and (X2A,X2B)
!
          IF (X1A.LT.X2B) THEN
            RLOW = MAX (X1A, X2A)
          ELSE
            GOTO 200
          ENDIF

          IF (X1B.GT.X2A) THEN
            RUPP = MIN (X1B, X2B)
          ELSE
            GOTO 200
          ENDIF

          IF (RUPP.LT.RLOW) THEN
            CALL MSGERR (3, 'interpolation error')
            WRITE (PRTEST, 140) I1, X1A, X1B, I2, X2A, X2B
 140        FORMAT (' I, XA, XB ', 2(I3, 2(1X,E12.4)))
          ELSE
            Y2(I2) = Y2(I2) + Y1(I1) * (RUPP-RLOW)
          ENDIF
 200    CONTINUE
!
!       Cell in X1 covers both ends of sector boundary
        IF (TWICE) THEN
          IF (X1A.LT.X2LO) THEN
             X1A = X1A + PERIOD
             X1B = X1B + PERIOD
          ENDIF
          IF (X1B.GT.X2HI) THEN
             X1A = X1A - PERIOD
             X1B = X1B - PERIOD
          ENDIF
          TWICE = .FALSE.
          GOTO 100
        ENDIF
 300  CONTINUE
!
      DO I2 = 1, N2
        IF (I2.EQ.1) THEN
          CELLEN = X2(2) - X2(1)
        ELSE IF (I2.EQ.N2) THEN
          CELLEN = X2(N2) - X2(N2-1)
        ELSE
          CELLEN = 0.5 * (X2(I2+1) - X2(I2-1))
        ENDIF
!       divide Y2 by cell length
        Y2(I2) = Y2(I2) / CELLEN
      ENDDO
      IF (ITEST.GE.160) THEN
        WRITE (PRTEST, 84) N1, N2
  84    FORMAT (' test CHGBAS ', 2I5)
        WRITE (PRTEST, 85) (X1(II), II = 1, N1)
        WRITE (PRTEST, 85) (Y1(II), II = 1, N1)
        WRITE (PRTEST, 85) (X2(II), II = 1, N2)
        WRITE (PRTEST, 85) (Y2(II), II = 1, N2)
  85    FORMAT (10 (1X,E10.3))
      ENDIF
!
      RETURN
      END
!
!********************************************************************
!                                                                   *
      REAL FUNCTION GAMMA(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)
!
      REAL XX, YY, ABIG                                                   40.00
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30./
      CALL STRACE (IENT, 'GAMMA')
      YY = GAMMLN(XX)
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMA = EXP(YY)
      RETURN
      END
!********************************************************************
!                                                                   *
      FUNCTION GAMMLN(XX)
!                                                                   *
!********************************************************************
!
!   Method:
!     function is copied from: Press et al., "Numerical Recipes"
!
      DOUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE WRSPEC (NREF, ACLOC)
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE OUTP_DATA                                                       40.13

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
!     40.00, 40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     new subroutine, update 40.00
!     40.03, Mar. 00: precision increased; 2 decimals more in output table
!     40.13, July 01: variable format using module OUTP_DATA
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Writing of action density spectrum in Swan standard format
!
!  3. METHOD
!
!
!  4. Argument variables
!
!       NREF    int    input    unit ref. number of output file
!       ACLOC   real   local    2-D spectrum or source term at one
!                               output location

      INTEGER, INTENT(IN) :: NREF
      REAL, INTENT(IN)    :: ACLOC(1:MDC,1:MSC)

!  5. Parameter variables
!
!  6. Local variables
!
!       ID      counter of spectral directions
!       IS      counter of spectral frequencies
!
      INTEGER :: ID, IS
!
!       EFAC    multiplication factor written to file
!
      REAL    :: EFAC
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       determine maximum value of ACLOC
!       if maximum = 0
!       then write 'ZERO' to file
!       else write 'FACTOR'
!            determine multiplication factor, write this to file
!            write values of ACLOC/factor to file
!       ----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0                                           40.13
      IF (LTRACE) CALL STRACE (IENT, 'WRSPEC')
!
!     first determine maximum energy density
      EFAC = 0.
      DO ID = 1, MDC
        DO IS = 1, MSC
          IF (ACLOC(ID,IS).GE.0.) THEN
            EFAC = MAX (EFAC, ACLOC(ID,IS))
          ELSE
            EFAC = MAX (EFAC, 10.*ABS(ACLOC(ID,IS)))
          ENDIF
        ENDDO
      ENDDO
      IF (EFAC .LE. 1.E-10) THEN
        WRITE (NREF, 12) 'ZERO'                                           40.00
  12    FORMAT (A4)
      ELSE
        EFAC = 1.01 * EFAC * 10.**(-DEC_SPEC)                             40.13
!       factor PI/180 introduced to account for change from rad to degr
!       factor 2*PI to account for transition from rad/s to Hz
        WRITE (NREF, 95) EFAC * 2. * PI**2 / 180.
  95    FORMAT ('FACTOR', /, E18.8)                                       40.13
        DO IS = 1, MSC
!         write spectral energy densities to file
          WRITE (NREF, FIX_SPEC) (NINT(ACLOC(ID,IS)/EFAC), ID=1,MDC)      40.13
        ENDDO
      ENDIF
      RETURN
!     end of subroutine WRSPEC
      END
!****************************************************************
!
      SUBROUTINE SWACC(AC2, AC2OLD, ACNRMS, ISSTOP, IDCMIN, IDCMAX)
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Sep. 02: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Determine some infinity norms meant for stop criterion
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     ACNRMS      array containing infinity norms
!     IDCMIN      integer array containing minimum counter of directions
!     IDCMAX      integer array containing maximum counter of directions
!     ISSTOP      maximum frequency counter in this sweep
!
      INTEGER IDCMIN(MSC), IDCMAX(MSC), ISSTOP
      REAL    AC2(MDC,MSC,MCGRD), AC2OLD(MDC,MSC), ACNRMS(2)
!
!  6. Local variables
!
!     DIFFAC:     difference between AC2 and AC2OLD
!     ID    :     counter of direction
!     IDDUM :     uncorrected counter of direction
!     IENT  :     number of entries
!     IS    :     counter of frequency
!
      INTEGER ID, IDDUM, IS, IENT
      REAL    DIFFAC
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWOMPU (in SWANCOM1)
!
! 12. Structure
!
!     determine infinity norms |ac2 - ac2old| and |ac2|
!     in selected sweep
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWACC')

      DO IS = 1, ISSTOP
         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
!           *** determine infinity norms |ac2 - ac2old| and |ac2|
!
            DIFFAC = ABS(AC2(ID,IS,KCGRD(1)) - AC2OLD(ID,IS))
            IF (DIFFAC.GT.ACNRMS(1)) ACNRMS(1) = DIFFAC
            IF (ABS(AC2(ID,IS,KCGRD(1))).GT.ACNRMS(2))
     &                              ACNRMS(2) = ABS(AC2(ID,IS,KCGRD(1)))

         END DO
      END DO

      RETURN
      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWTSTA (ITIMER)
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE TIMECOMM                                                        40.41
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG!
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Start timing
!TIMG!
!TIMG!  3. Method
!TIMG!
!TIMG!     Get cpu and wall-clock times and store
!TIMG!
!TIMG!  4. Argument variables
!TIMG!
!TIMG!     ITIMER      number of timer to be used
!TIMG!
!TIMG      INTEGER :: ITIMER
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     C     :     clock count of processor
!TIMG!     I     :     index in LISTTM, loop variable
!TIMG!     IFOUND:     index in LISTTM, location of ITIMER
!TIMG!     IFREE :     index in LISTTM, first free position
!TIMG!     M     :     maximum clock count
!TIMG!     R     :     number of clock counts per second
!TIMG!F95!     TIMER :     current real cpu-time
!TIMG!     TIMER1:     current cpu-time used
!TIMG!     TIMER2:     current wall-clock time used
!TIMG!
!TIMG      INTEGER          :: I, IFOUND, IFREE
!TIMG      INTEGER          :: C, R, M
!TIMG!F95      REAL             :: TIMER
!TIMG      DOUBLE PRECISION :: TIMER1, TIMER2
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!F95!     CPU_TIME         Returns real value from cpu-time clock
!TIMG!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN, SWCOMP, SWOMPU
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Get and store the cpu and wall-clock times
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG
!TIMG!
!TIMG!     --- check whether a valid timer number is given
!TIMG!
!TIMG      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
!TIMG         WRITE(PRINTF,*) 'SWTSTA: ITIMER out of range: ',
!TIMG     &                   ITIMER, 1, NSECTM
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- check whether timing for ITIMER was started already,
!TIMG!         also determine first free location in LISTTM
!TIMG!
!TIMG      IFOUND=0
!TIMG      IFREE =0
!TIMG      I     =0
!TIMG 100  IF (I.LT.LASTTM .AND. (IFOUND.EQ.0 .OR. IFREE.EQ.0)) THEN
!TIMG         I=I+1
!TIMG         IF (LISTTM(I).EQ.ITIMER) THEN
!TIMG            IFOUND=I
!TIMG         END IF
!TIMG         IF (IFREE.EQ.0 .AND. LISTTM(I).EQ.-1) THEN
!TIMG            IFREE =I
!TIMG         END IF
!TIMG         GOTO 100
!TIMG      END IF
!TIMG
!TIMG      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0 .AND. LASTTM.LT.MXTIMR) THEN
!TIMG         LASTTM=LASTTM+1
!TIMG         IFREE =LASTTM
!TIMG      END IF
!TIMG!
!TIMG!     --- produce warning if found in the list
!TIMG!
!TIMG      IF (IFOUND.GT.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTA: warning: previous timing for section ',
!TIMG     &      ITIMER,' not closed properly/will be ignored.'
!TIMG      END IF
!TIMG!
!TIMG!     --- produce error if not found and no free position available
!TIMG!
!TIMG      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTA: maximum number of simultaneous timers',
!TIMG     &      ' exceeded:',MXTIMR
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- register ITIMER in appropriate location of LISTTM
!TIMG!
!TIMG      IF (IFOUND.EQ.0) THEN
!TIMG         IFOUND=IFREE
!TIMG      END IF
!TIMG      LISTTM(IFOUND)=ITIMER
!TIMG!
!TIMG!     --- get current cpu/wall-clock time and store in TIMERS
!TIMG!
!TIMG      TIMER1=0D0
!TIMG!F95      CALL CPU_TIME (TIMER)
!TIMG!F95      TIMER1=DBLE(TIMER)
!TIMG      CALL SYSTEM_CLOCK (C,R,M)
!TIMG      TIMER2=DBLE(C)/DBLE(R)
!TIMG
!TIMG      TIMERS(IFOUND,1)=TIMER1
!TIMG      TIMERS(IFOUND,2)=TIMER2
!TIMG
!TIMG      RETURN
!TIMG      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWTSTO (ITIMER)
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE TIMECOMM                                                        40.41
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG!
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Stop timing
!TIMG!
!TIMG!  3. Method
!TIMG!
!TIMG!     Get cpu and wall-clock times and store
!TIMG!
!TIMG!  4. Argument variables
!TIMG!
!TIMG!     ITIMER      number of timer to be used
!TIMG!
!TIMG      INTEGER :: ITIMER
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     C     :     clock count of processor
!TIMG!     I     :     index in LISTTM, loop variable
!TIMG!     IFOUND:     index in LISTTM, location of ITIMER
!TIMG!     M     :     maximum clock count
!TIMG!     R     :     number of clock counts per second
!TIMG!F95!     TIMER :     current real cpu-time
!TIMG!     TIMER1:     current cpu-time used
!TIMG!     TIMER2:     current wall-clock time used
!TIMG!
!TIMG      INTEGER          :: I, IFOUND
!TIMG      INTEGER          :: C, R, M
!TIMG!F95      REAL             :: TIMER
!TIMG      DOUBLE PRECISION :: TIMER1, TIMER2
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!F95!     CPU_TIME         Returns real value from cpu-time clock
!TIMG!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN, SWCOMP, SWOMPU
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Get and store the cpu and wall-clock times
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG
!TIMG!
!TIMG!     --- check whether a valid timer number is given
!TIMG!
!TIMG      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
!TIMG         WRITE(PRINTF,*) 'SWTSTO: ITIMER out of range: ',
!TIMG     &                   ITIMER, 1, NSECTM
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- check whether timing for ITIMER was started already,
!TIMG!         also determine first free location in LISTTM
!TIMG!
!TIMG      IFOUND=0
!TIMG      I     =0
!TIMG 100  IF (I.LT.LASTTM .AND. IFOUND.EQ.0) THEN
!TIMG         I=I+1
!TIMG         IF (LISTTM(I).EQ.ITIMER) THEN
!TIMG            IFOUND=I
!TIMG         END IF
!TIMG         GOTO 100
!TIMG      END IF
!TIMG!
!TIMG!     --- produce error if not found
!TIMG!
!TIMG      IF (IFOUND.EQ.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTO: section ',ITIMER,' not found',
!TIMG     &      ' in list of active timings'
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- get current cpu/wall-clock time
!TIMG!
!TIMG      TIMER1=0D0
!TIMG!F95      CALL CPU_TIME (TIMER)
!TIMG!F95      TIMER1=DBLE(TIMER)
!TIMG      CALL SYSTEM_CLOCK (C,R,M)
!TIMG      TIMER2=DBLE(C)/DBLE(R)
!TIMG!
!TIMG!     --- calculate elapsed time since start of timing,
!TIMG!         store in appropriate location in DCUMTM,
!TIMG!         increment number of timings for current section
!TIMG!
!TIMG      DCUMTM(ITIMER,1)=DCUMTM(ITIMER,1)+(TIMER1-TIMERS(IFOUND,1))
!TIMG      DCUMTM(ITIMER,2)=DCUMTM(ITIMER,2)+(TIMER2-TIMERS(IFOUND,2))
!TIMG      NCUMTM(ITIMER)  =NCUMTM(ITIMER)+1
!TIMG!
!TIMG!     --- free appropriate location of LISTTM,
!TIMG!         adjust last occupied position of LISTTM
!TIMG!
!TIMG      IF (IFOUND.GT.0) THEN
!TIMG         LISTTM(IFOUND)=-1
!TIMG      END IF
!TIMG 200  IF (LASTTM.GT.1 .AND. LISTTM(LASTTM).EQ.-1) THEN
!TIMG         LASTTM=LASTTM-1
!TIMG         GOTO 200
!TIMG      END IF
!TIMG      IF (LISTTM(LASTTM).EQ.-1) LASTTM=0
!TIMG
!TIMG      RETURN
!TIMG      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWPRTI
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG      USE TIMECOMM                                                        40.41
!TIMG      USE M_PARALL                                                        40.31
!TIMG!PUN      USE SIZES, ONLY: MNPROC, MYPROC                                     40.95
!TIMG
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.30: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Print timings info
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     IDEBUG:     level of timing output requested:
!TIMG!                 0 - no output for detailed timings
!TIMG!                 1 - aggregate output for detailed timings
!TIMG!                 2 - complete output for all detailed timings
!TIMG!     IENT  :     number of entries
!TIMG!     J     :     loop counter
!TIMG!     K     :     loop counter
!TIMG!     MYPRC :     own process number
!TIMG!     TABLE :     array for computing aggregate cpu- and wallclock-times
!TIMG!
!TIMG      INTEGER          :: IENT, J, K, IDEBUG, MYPRC
!TIMG      DOUBLE PRECISION :: TABLE(30,2)
!TIMG      PARAMETER (IDEBUG=0)
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!     STRACE           Tracing routine for debugging
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN (in SWANMAIN)
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Compile table with overview of cpu/wall clock time used in
!TIMG!     important parts of SWAN and write to PRINT file
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG      SAVE IENT
!TIMG      DATA IENT/0/
!TIMG      IF (LTRACE) CALL STRACE (IENT,'SWPRTI')
!TIMG!
!TIMG      MYPRC = INODE
!TIMG!PUN      IF ( MNPROC>1 ) MYPRC = MYPROC
!TIMG!
!TIMG!     --- compile table with overview of cpu/wall clock time used in
!TIMG!         important parts of SWAN and write to PRINT file
!TIMG!
!TIMG      IF ( ITEST.GE.1 .OR. IDEBUG.GE.1 ) THEN
!TIMG!
!TIMG!        --- initialise table to zero
!TIMG!
!TIMG         DO K = 1, 30
!TIMG            DO J = 1, 2
!TIMG               TABLE(K,J) = 0D0
!TIMG            END DO
!TIMG         END DO
!TIMG!
!TIMG!        --- compute times for basic blocks
!TIMG!
!TIMG         DO J = 1, 2
!TIMG!
!TIMG!           --- total run-time
!TIMG!
!TIMG            TABLE(1,J) = DCUMTM(1,J)
!TIMG!
!TIMG!           --- initialisation, reading, preparation:
!TIMG!
!TIMG            DO K = 2, 7
!TIMG               TABLE(2,J) = TABLE(2,J) + DCUMTM(K,J)
!TIMG            END DO
!TIMG!
!TIMG!           --- domain decomposition:
!TIMG!
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(211,J)
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(212,J)
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(201,J)
!TIMG!
!TIMG!           --- total calculation including communication:
!TIMG!
!TIMG            TABLE(3,J) = TABLE(3,J) + DCUMTM(8,J)
!TIMG!
!TIMG!           --- output:
!TIMG!
!TIMG            TABLE(5,J) = TABLE(5,J) + DCUMTM(9,J)
!TIMG!
!TIMG!           --- exchanging data:
!TIMG!
!TIMG            TABLE(7,J) = TABLE(7,J) + DCUMTM(213,J)
!TIMG!
!TIMG!           --- solving system:
!TIMG!
!TIMG            TABLE(9,J) = TABLE(9,J) + DCUMTM(119,J)
!TIMG            TABLE(9,J) = TABLE(9,J) + DCUMTM(120,J)
!TIMG!
!TIMG!           --- global reductions:
!TIMG!
!TIMG            TABLE(10,J) = TABLE(10,J) + DCUMTM(202,J)
!TIMG!
!TIMG!           --- collecting data:
!TIMG!
!TIMG            TABLE(11,J) = TABLE(11,J) + DCUMTM(214,J)
!TIMG!
!TIMG!           --- setup:
!TIMG!
!TIMG            TABLE(12,J) = TABLE(12,J) + DCUMTM(106,J)
!TIMG!
!TIMG!           --- propagation velocities:
!TIMG!
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(111,J)
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(113,J)
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(114,J)
!TIMG!
!TIMG!           --- x-y advection:
!TIMG!
!TIMG            TABLE(15,J) = TABLE(15,J) + DCUMTM(140,J)
!TIMG!
!TIMG!           --- sigma advection:
!TIMG!
!TIMG            TABLE(16,J) = TABLE(16,J) + DCUMTM(141,J)
!TIMG!
!TIMG!           --- theta advection:
!TIMG!
!TIMG            TABLE(17,J) = TABLE(17,J) + DCUMTM(142,J)
!TIMG!
!TIMG!           --- wind:
!TIMG!
!TIMG            TABLE(18,J) = TABLE(18,J) + DCUMTM(132,J)
!TIMG!
!TIMG!           --- whitecapping:
!TIMG!
!TIMG            TABLE(19,J) = TABLE(19,J) + DCUMTM(133,J)
!TIMG!
!TIMG!           --- bottom friction:
!TIMG!
!TIMG            TABLE(20,J) = TABLE(20,J) + DCUMTM(130,J)
!TIMG!
!TIMG!           --- wave breaking:
!TIMG!
!TIMG            TABLE(21,J) = TABLE(21,J) + DCUMTM(131,J)
!TIMG!
!TIMG!           --- quadruplets:
!TIMG!
!TIMG            TABLE(22,J) = TABLE(22,J) + DCUMTM(135,J)
!TIMG!
!TIMG!           --- triads:
!TIMG!
!TIMG            TABLE(23,J) = TABLE(23,J) + DCUMTM(134,J)
!TIMG!
!TIMG!           --- limiter:
!TIMG!
!TIMG            TABLE(24,J) = TABLE(24,J) + DCUMTM(122,J)
!TIMG!
!TIMG!           --- rescaling:
!TIMG!
!TIMG            TABLE(25,J) = TABLE(25,J) + DCUMTM(121,J)
!TIMG!
!TIMG!           --- reflections:
!TIMG!
!TIMG            TABLE(26,J) = TABLE(26,J) + DCUMTM(136,J)
!TIMG!
!TIMG!           --- diffraction:
!TIMG!
!TIMG            TABLE(27,J) = TABLE(27,J) + DCUMTM(137,J)
!TIMG!
!TIMG!           --- vegetation:
!TIMG!
!TIMG            TABLE(29,J) = TABLE(29,J) + DCUMTM(139,J)
!TIMG
!TIMG         END DO
!TIMG!
!TIMG!        --- add up times for some basic blocks
!TIMG!
!TIMG         DO J = 1, 2
!TIMG!
!TIMG!           --- total calculation:
!TIMG!
!TIMG            TABLE(3,J) = TABLE(3,J) - TABLE( 7,J)
!TIMG            TABLE(3,J) = TABLE(3,J) - TABLE(10,J)
!TIMG            IF ( TABLE(3,J).LT.0D0 ) TABLE(3,J) = 0D0
!TIMG!
!TIMG!           --- total communication:
!TIMG!                * exchanging data
!TIMG!                * global reductions
!TIMG!                * collecting data
!TIMG!
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE( 7,J)
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE(10,J)
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE(11,J)
!TIMG!
!TIMG!           --- total propagation:
!TIMG!                * velocities and derivatives
!TIMG!
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(14,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(15,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(16,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(17,J)
!TIMG!
!TIMG!           --- sources:
!TIMG!                * wind, whitecapping, friction, breaking,
!TIMG!                * quadruplets, triads, limiter, rescaling,
!TIMG!                * reflections
!TIMG!
!TIMG            DO K = 18, 26
!TIMG               TABLE(8,J) = TABLE(8,J) + TABLE(K,J)
!TIMG            END DO
!TIMG!
!TIMG!                * diffraction
!TIMG!
!TIMG            TABLE(8,J) = TABLE(8,J) + TABLE(27,J)
!TIMG!
!TIMG!                * vegetation
!TIMG!
!TIMG            TABLE(8,J) = TABLE(8,J) + TABLE(29,J)
!TIMG!
!TIMG!           --- other computing:
!TIMG!
!TIMG            TABLE(13,J) = TABLE(13,J) + TABLE( 3,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 6,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 8,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 9,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE(12,J)
!TIMG            IF ( TABLE(13,J).LT.0D0 ) TABLE(13,J) = 0D0
!TIMG
!TIMG         END DO
!TIMG!
!TIMG!        --- print CPU-times used in important parts of SWAN
!TIMG!
!TIMG         WRITE(PRINTF,'(/)')
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,111) MYPRC
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,112) MYPRC
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,115) MYPRC,'total time:'       ,(TABLE(1,J),J=1,2)
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,115) MYPRC,'total pre-processing:',
!TIMG     &                                                (TABLE(2,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'total calculation:',(TABLE(3,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'total communication:',
!TIMG     &                                                (TABLE(4,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'total post-processing:',
!TIMG     &                                                (TABLE(5,j),j=1,2)
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,113) MYPRC
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,115) MYPRC,'calc. propagation:',(TABLE(6,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'exchanging data:'  ,(TABLE(7,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'calc. sources:'    ,(TABLE(8,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'solving system:'   ,(TABLE(9,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'reductions:'      ,(TABLE(10,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'collecting data:' ,(TABLE(11,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'calc. setup:'     ,(TABLE(12,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'other computing:' ,(TABLE(13,j),j=1,2)
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,114) MYPRC
!TIMG         WRITE(PRINTF,110) MYPRC
!TIMG         WRITE(PRINTF,115) MYPRC,'prop. velocities:',(TABLE(14,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'x-y advection:'   ,(TABLE(15,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'sigma advection:' ,(TABLE(16,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'theta advection:' ,(TABLE(17,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'wind:'            ,(TABLE(18,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'whitecapping:'    ,(TABLE(19,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'bottom friction:' ,(TABLE(20,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'vegetation:'      ,(TABLE(29,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'wave breaking:'   ,(TABLE(21,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'quadruplets:'     ,(TABLE(22,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'triads:'          ,(TABLE(23,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'limiter:'         ,(TABLE(24,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'rescaling:'       ,(TABLE(25,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'reflections:'     ,(TABLE(26,j),j=1,2)
!TIMG         WRITE(PRINTF,115) MYPRC,'diffraction:'     ,(TABLE(27,j),j=1,2)
!TIMG
!TIMG      END IF
!TIMG
!TIMG      IF ( IDEBUG.GE.2 ) THEN
!TIMG         WRITE(PRINTF,120) MYPRC
!TIMG         DO J = 1, NSECTM
!TIMG            IF (NCUMTM(J).GT.0)
!TIMG     &         WRITE(PRINTF,121) MYPRC,J,DCUMTM(J,1),DCUMTM(J,2),
!TIMG     &                           NCUMTM(J)
!TIMG         END DO
!TIMG      END IF
!TIMG
!TIMG  110 FORMAT(i3,1x,'#')
!TIMG  111 FORMAT(i3,' # Details on timings of the simulation:')
!TIMG  112 FORMAT(i3,1x,'#',26x,'cpu-time',1x,'wall-clock')
!TIMG  113 FORMAT(i3,' # Splitting up calc. + comm. times:')
!TIMG  114 FORMAT(i3,' # Overview source contributions:')
!TIMG  115 FORMAT(i3,1x,'#',1x,a22,2f11.2)
!TIMG
!TIMG  120 FORMAT(/,i3,' #    item     cpu-time    real time     count')
!TIMG  121 FORMAT(i3,1x,'#',4x,i4,2f13.4,i10)
!TIMG
!TIMG      RETURN
!TIMG      END
!****************************************************************
!
      SUBROUTINE TXPBLA(TEXT,IF,IL)
!
!****************************************************************
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     determines the position of the first and the last non-blank
!     (or non-tabulator) character in the text-string
!
!  4. Argument variables
!
!     IF          position of the first non-blank character in TEXT
!     IL          position of the last non-blank character in TEXT
!     TEXT        text string
!
      INTEGER IF, IL
      CHARACTER*(*) TEXT
!
!  6. Local variables
!
!     FOUND :     TEXT is found or not
!     ITABVL:     integer value of tabulator character
!     LENTXT:     length of TEXT
!
      INTEGER LENTXT, ITABVL
      LOGICAL FOUND
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
!DOS      ITABVL = ICHAR('	')
      ITABVL = ICHAR('	')
      LENTXT = LEN (TEXT)
      IF = 1
      FOUND = .FALSE.
  100 IF (IF .LE. LENTXT .AND. .NOT. FOUND) THEN
         IF (.NOT. (TEXT(IF:IF) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IF:IF)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ELSE
            IF = IF + 1
         ENDIF
         GOTO 100
      ENDIF
      IL = LENTXT + 1
      FOUND = .FALSE.
  200 IF (IL .GT. 1 .AND. .NOT. FOUND) THEN
         IL = IL - 1
         IF (.NOT. (TEXT(IL:IL) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IL:IL)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ENDIF
         GOTO 200
      ENDIF

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION INTSTR ( IVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer to string
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!
      INTEGER IVAL
!
!  6. Local variables
!
!     CVAL  :     character represented an integer of mantisse
!     I     :     counter
!     IPOS  :     position in mantisse
!     IQUO  :     whole quotient
!
      INTEGER I, IPOS, IQUO
      CHARACTER*1, ALLOCATABLE :: CVAL(:)
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IPOS = 1
 100  CONTINUE
      IF (IVAL/10**IPOS.GE.1.) THEN
         IPOS = IPOS + 1
         GO TO 100
      END IF
      ALLOCATE(CVAL(IPOS))

      DO I=IPOS,1,-1
         IQUO=IVAL/10**(I-1)
         CVAL(IPOS-I+1)=CHAR(INT(IQUO)+48)
         IVAL=IVAL-IQUO*10**(I-1)
      END DO

      WRITE (INTSTR,*) (CVAL(I), I=1,IPOS)

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION NUMSTR ( IVAL, RVAL, FORM )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Convert integer or real to string with given format
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!     FORM        given format
!     RVAL        real to be converted
!
      INTEGER   IVAL
      REAL      RVAL
      CHARACTER FORM*20
!
!  6. Local variables
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IF ( IVAL.NE.INAN ) THEN
         WRITE (NUMSTR,FORM) IVAL
      ELSE IF ( RVAL.NE.RNAN ) THEN
         WRITE (NUMSTR,FORM) RVAL
      ELSE
         NUMSTR = ''
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPI ( IARR1, IARR2, LENGTH )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies integer array IARR1 to IARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IARR1       source array
!     IARR2       target array
!     LENGTH      array length
!
      INTEGER LENGTH
      INTEGER IARR1(LENGTH), IARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
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
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPI')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array IARR1 to IARR2

      DO 10 I = 1, LENGTH
         IARR2(I) = IARR1(I)
  10  CONTINUE

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPR ( ARR1, ARR2, LENGTH )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies real array ARR1 to ARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     ARR1        source array
!     ARR2        target array
!     LENGTH      array length
!
      INTEGER LENGTH
      REAL    ARR1(LENGTH), ARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
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
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPR')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array ARR1 to ARR2

      DO 10 I = 1, LENGTH
         ARR2(I) = ARR1(I)
  10  CONTINUE

      RETURN
      END
