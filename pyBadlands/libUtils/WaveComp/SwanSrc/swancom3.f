!
!     SWAN/COMPU   file 3 of 5
!
!     PROGRAM SWANCOM3.FOR
!
!     This file SWANCOM3 of the main program SWAN program
!     includes the next subroutines (mainly subroutines for
!     the source terms for generation of wave energy ) :
!
!     WNDPAR  (DOLPHIN-B formulations for the SWAN model for a
!             first- or a second generation first guess of the spectrum)
!     WINDP1 (computation of variables derived from the wind such as
!             mean wind velocity, mean wind direction, minimum counter
!             for the wind, maximum counter for the wind, wind friction
!             velocity and the Pierson Moskowitz frequency )
!     WINDP2 (computation of wind sea energy spectrum necessary for
!             the second generation wind growth model)
!     WINDP3 (limit the energy spectrum in the case of a first or
!             second generation wind growth model)
!     SWIND0 (linear input term Cavaleri and Malanotte Rizolli (1981)
!     SWIND3 (third generation wind growth model (Snyder et al. 1981;
!             Komen et al., 1984)
!     SWIND4 (third generation wind growth model (Janssen, 1989,1991)
!     SWIND5 (third generation wind growth model according to the
!             expression of Yan (1987) (especially derived for a
!             frequency range that extend to the high frequencies
!
!****************************************************************
!
      SUBROUTINE WNDPAR (ISSTOP,IDWMIN,IDWMAX,IDCMIN,IDCMAX,              32.06
     &                   DEP2  ,WIND10,GENC0 ,GENC1 ,                     40.85 32.06
     &                   THETAW,AC2   ,KWAVE ,IMATRA,IMATDA,              32.06
     &                   SPCSIG,CGO   ,ALIMW ,GROWW ,ETOTW ,              32.06
     &                   PLWNDA,PLWNDB,SPCDIR,ITER            )           32.06
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE                                                       30.82
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
!     30.75: IJsbrand Haagsma (bug fix)
!     30.82: IJsbrand Haagsma
!     32.06: Roeland Ris
!     40.00: Nico Booij (Nonstationary boundary conditions)
!     40.41: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!            Jan. 97: New subroutine (Roeland Ris)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: argument list of WINDP2 changed
!     30.75, Mar. 98: Set FPM=SIGPKD, due to change in argument list of WINDP2
!     40.00, July 98: argument list of WINDP2 changed
!     30.82, Oct. 98: Updated description of several variables
!     30.82, Apr. 99: Dimensioning KCGRD corrected
!     32.06, June 99: Reformulated directional spreading for first guess
!     30.82, June 99: Implicit none added; all variables declared
!     40.41, Aug. 04: COS(THETA-THETAW) replaced by sumrule to make it cheaper
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.85, Aug. 08: store wind input for output purposes
!
!  2. Purpose
!
!     Computation of the wind input source term with formulations
!     of a first-generation model (constant porportionality coefficient)
!     and a second-generation model (proportionality coefficient depends
!     on the energy in the wind sea part of the spectrum). The
!     expressions are from Holthuijsen and de Boer (1988) and from
!     the DOLPHIN-B model (Holthuijsen and Booij). During the
!     implementation of the terms modifications to the code have been
!     made after personal communications with Holthuijsen and Booij.
!
!  3. Method
!
!     The source term of the following nature:
!
!     S = A + B E          for E < Elim   | t - tw | < pi/2
!
!         (Elim-E)
!     S = --------         for E > Elim   | t - tw | < pi/2
!            TAU
!
!     S = 0                for E > Elim   | t - tw | > pi/2
!
!     in which the terms A and B are:
!
!         [cf10]          2        2        2              4
!     A = ------ pi (1./g)  [rhoaw]  [cdrag]  (U cos(t-tw))
!          2 pi
!
!     and:
!                              U cos(t-tw)              s
!     B = 5 [cf20] [rhoaw] f { ----------- -  [cf30] } ----
!                                  Cph                 2 pi
!
!     The coefficient TAU in the relaxation model is given by:
!                           2
!                   / 2 pi \      g
!     TAU = [cf40] | -----  | ------------
!                   \  s    /  U cos(d-dw)
!
!     The limiting spectrum is given by:
!
!                    -3
!             ALPHA K                 s   -4    2     2
!     Elim = ----------- exp ( -5/4 (----)   ) --- cos ( t - tw )
!             2  Cg                  spm        pi
!
!     in which:
!
!        ALPHA   wind sea and/or depth dependent proportionality
!                coefficient which controls the energy scale of the
!                limiting spectrum.
!              * In the first-generation model ALPHA is a constant
!                equal to 0.0081 (fully developed)
!              * In the second-generation model ALPHA depends on the
!                energy in the wind sea part of the spectrum. ALPHA
!                is calculated here by:
!                                                [cf60]
!                            ALPHA = [cf50] * Edml
!
!        spm     adapted Pierson-Moskowitz (1964) peak frequency
!
!     The total non-dimensional energy in the wind sea part of the
!     spectrum is calculated by (see subroutine WINDP2):
!
!                   2
!               grav  * ETOTW
!       Edml =  -------------
!                       4
!                 wind10
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
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!  6. Local variables
!
!     IENT  : Number of entries into this subroutine
!
      INTEGER IENT
!
      REAL  :: FPM    ! Pierson-Moskowitz frequency
      REAL  :: SWIND_EXP, SWIND_IMP    ! explicit and implicit part of wind source
!
!        INTEGERS:
!        ---------
!        IDWMIN      Minimum counter for spectral wind direction
!        IDWMAX      Maximum counter for spectral wind direction
!        IX          Counter of gridpoint in x-direction
!        IY          Counter of gridpoint in y-direction
!        IS          Counter of frequency bin
!        ISSTOP      Countrer for the maximum frequency of all directions
!        IDDUM       Dummy counter
!        ID          Counter of directional distribution
!        IDWMIN/IDWMAX  Minimum / maximum counter in wind sector (180 degrees)
!
!        REALS:
!        ---------
!        ALPM        Coefficient for overshoot at deep water
!        ALPMD       Coefficient for overshoot corrected for shallow
!                    water using expression of Bretschnieder (1973)
!        ALIMW       limiting spectrum in terms of action density
!        ARG1, ARG2  Exponent
!        CDRAG       Wind drag coefficient
!        DND         Nondimensional depth
!        DTHETA      Difference in rad between wave and wind direction
!        EDML        Dimensionless energy
!        ETOTW       Total energy of the wind sea part of the spectrum
!        RHOAW       Density of air devided by the density of water
!        SIGPK       Peak frequency in terms of rad /s
!        SIGPKD      Adapted peak frequency for shallow water
!        TAU         Variable for the wind growth equation
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        TWOPI       Two times pi
!        WIND10      Velocity of the relative wind vector
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!        ALIMW     2D    Limiting action density spectrum
!        DEP2      1D    Depth
!        CGO       2D    Group velocity
!        KWAVE     2D    Wave number
!        LOGSIG    1D    Logaritmic distribution of frequency
!        IMATRA    2D    Coefficients of right hand side of vector
!        IMATDA    2D    Coefficients of the diagonal
!        PLWNDA    3D    Values of source term for test point
!        PLWNDB    3D    Values of source term for test point
!        SPCDIR    1D    Spectral direction of wave component
!        IDCMIN    1D    Minimum counter
!        IDCMAX    1D    Maximum counter in directional space
!        GROWW     2D    Aux. array to determine whether there are
!                        wave generation conditions
!
!        PWIND(1)  = CF10     188.0
!        PWIND(2)  = CF20     0.59
!        PWIND(3)  = CF30     0.12
!        PWIND(4)  = CF40     250.0
!        PWIND(5)  = CF50     0.0023
!        PWIND(6)  = CF60    -0.2233
!        PWIND(7)  = CF70     0.       (not used)
!        PWIND(8)  = CF80    -0.56     (not used)
!        PWIND(9)  = RHOAW    0.00125  (density air / density water)
!        PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!        PWIND(11) = CDRAG    0.0012   (drag coefficient)
!        PWIND(12) = UMIN     1.0      (minimum wind velocity)
!        PWIND(13) = PMLM     0.13     (  )
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
!
!     6. SUBROUTINES USED
!
!        WINDP2     (compute the total energy in the wind sea part of
!                     the spectrum). Subroutine WINDP2 is called in
!                     SWANCOM1 in subroutine SOURCE
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!   --------------------------------------------------------------
!   Calculate the adapted peak frequency
!   --------------------------------------------------------------
!   If first-generation model
!     alpha is constant
!   else
!     Calculate energy in wind sea part of spectrum ETOTW
!     Calculate alpha on the basis of ETOTW
!   end
!   --------------------------------------------------------------
!   Take depth effects into account for alpha
!   --------------------------------------------------------------
!   For each frequency and direction
!     compute limiting spectrum and determine whether there is
!     grow or decay
!   end
!   --------------------------------------------------------------
!   Do for each frequency and direction
!     If wind-wave generation conditions are present
!       calculate A + B E
!     else if energy is larger than limiting spectrum
!       calculate dissipation rate with relaxation model
!     endif
!   enddo
!   --------------------------------------------------------------
!   Store results in matrix (IMATRA or IMATDA)
!   --------------------------------------------------------------
!   End of the subroutine WNDPAR
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER  IS    ,ID    ,ITER  ,
     &         IDWMIN,IDWMAX,IDDUM ,ISSTOP
!
      REAL     WIND10,THETA ,THETAW,EDML  ,ARG1  ,ARG2  ,
     &         ALPM  ,ALPMD ,TEMP1 ,TEMP2 ,FACTA ,FACTB ,
     &         ADUM  ,BDUM  ,CINV  ,SIGTPI,SIGMA ,TWOPI ,TAUINV,
     &         SIGPK ,SIGPKD,DND   ,ETOTW ,ALIM1D,
     &         CTW   ,STW   ,COSDIF,                                      40.41
     &         DIRDIS,AC2CEN,DTHETA
!
      REAL  :: AC2(MDC,MSC,MCGRD)
      REAL  :: ALIMW(MDC,MSC)
      REAL  :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  :: KWAVE(MSC,MICMAX)                                          40.22
      REAL  :: PLWNDA(MDC,MSC,NPTST)                                      40.00
      REAL  :: PLWNDB(MDC,MSC,NPTST)
      REAL  :: GENC0(MDC,MSC,MGENR)                                       40.85
      REAL  :: GENC1(MDC,MSC,MGENR)                                       40.85
      REAL  :: DEP2(MCGRD)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  :: CGO(MSC,MICMAX)                                            40.22
!
      INTEGER  IDCMIN(MSC)           ,
     &         IDCMAX(MSC)
!
      LOGICAL  GROWW(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'WNDPAR')
!
!     *** initialization of arrays ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          GROWW(ID,IS) = .FALSE.
          ALIMW(ID,IS) = 0.
        ENDDO
      ENDDO
!
!     *** calculate the adapted shallow water peak frequency         ***
!     *** according to Bretschneider (1973) using the nondimensional ***
!     *** depth DND                                                  ***
!
      TWOPI  = 2. * PI
      DND    = MIN( 50. , GRAV * DEP2(KCGRD(1)) / WIND10**2 )
      SIGPK  = TWOPI * 0.13 * GRAV / WIND10
      SIGPKD = SIGPK / TANH(0.833*DND**0.375)
      FPM    = SIGPKD                                                     30.75
      CTW    = COS(THETAW)                                                40.41
      STW    = SIN(THETAW)                                                40.41
!
      IF ( IWIND .EQ. 1 ) THEN
!
!       *** first generation model ***
!
        ALPM = 0.0081
!
      ELSE IF (IWIND .EQ. 2 ) THEN
!
!       *** second generation model ***
!
!       *** Determine the proportionality constant alpha on the basis ***
!       *** of the total energy in the wind sea part of the spectrum  ***
!       *** output of subroutine (WINDP2) is ETOTW                    ***
!
        CALL WINDP2 (IDWMIN  ,IDWMAX  ,SIGPKD  ,FPM     ,
     &               ETOTW   ,
     &               AC2     ,SPCSIG  ,         WIND10               )    40.00

        EDML = MIN ( PWIND(10) , (GRAV**2 * ETOTW) / WIND10**4 )
        EDML = MAX ( 1.E-25 , EDML )
!
        ARG1 = ABS(PWIND(6))
        ALPM = MAX( 0.0081, (PWIND(5) * (1./EDML)**ARG1) )
!
      ENDIF
!
!     *** Take into account depth effects for proportionality ***
!     *** constant alpha through the nondimensional depth DND ***
!
      ALPMD  = 0.0081 + ( 0.013 - 0.0081 ) * EXP ( -1. * DND )
      ALPM   = MIN ( 0.155  ,  MAX ( ALPMD , ALPM ) )
!
!     *** Calculate the limiting spectrum in terms of action density   ***
!     *** for the wind sea part (centered around the local wind        ***
!     *** direction). For conversion of f^-5 --> k^-3 and coefficients ***
!     *** see Kitaigorodskii et al. 1975                               ***
!
      DO IS = 1, ISSTOP
        TEMP1  = ALPM / ( 2. * KWAVE(IS,1)**3 * CGO(IS,1) )
        ARG2   = MIN ( 2. , SIGPKD / SPCSIG(IS) )                         30.72
        TEMP2  = EXP ( (-5./4.) * ARG2**4 )
        ALIM1D = TEMP1 * TEMP2 / SPCSIG(IS)                               30.72
        DO IDDUM = IDWMIN, IDWMAX
          ID     = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          THETA  = SPCDIR(ID,1)                                           30.82
          COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                    40.41
!
!     For better convergence the first guess of the directional spreading 32.06
!     is modified in third generation mode. The new formulation better    32.06
!     fits the directional spreading of the deep water growth curves.     32.06
!
          IF ((ITER.EQ.1).AND.(IGEN.EQ.3)) THEN                           32.06
            DIRDIS = 0.434917 * (MAX(0., COSDIF))**0.6                    40.41 32.06
          ELSE                                                            32.06
            DIRDIS = (2./PI) * COSDIF**2                                  40.41
          END IF                                                          32.06
!
          ALIMW(ID,IS) = ALIM1D * DIRDIS
          AC2CEN       = AC2(ID,IS,KCGRD(1))
          IF ( AC2CEN .LE. ALIMW(ID,IS) ) THEN
            GROWW(ID,IS) = .TRUE.
          ELSE
            GROWW(ID,IS) = .FALSE.
          ENDIF
        ENDDO
!       *** test output ***
        IF ( TESTFL .AND. ITEST .GE. 10 ) THEN
          WRITE(PRINTF,2002) IS, SPCSIG(IS), KWAVE(IS,1), CGO(IS,1)       30.72
 2002     FORMAT(' WNDPAR: IS SPCSIG KWAVE CGO :',I3,3E12.4)              30.72
          WRITE(PRINTF,2003) TEMP1, TEMP2, ARG2
 2003     FORMAT(' WNDPAR: TEMP1 TEMP2 ARG2    :',3X,3E12.4)
        END IF
      ENDDO
!
!     *** Calculate the wind input (linear term A and exponential  ***
!     *** term B) in wave generating conditions or disspation term ***
!     *** if energy in bin is larger than limiting spectrum        ***
!
!
      FACTA = PWIND(1) * PI * PWIND(9)**2 * PWIND(11)**2 / GRAV**2
!
      DO IS = 1, ISSTOP
        SIGMA   = SPCSIG(IS)                                              30.72
        SIGTPI  = SIGMA * TWOPI
        CINV    = KWAVE(IS,1) / SIGMA
        FACTB = PWIND(2) * PWIND(9) * SIGMA / TWOPI                       34.00
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID     = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          DTHETA = SPCDIR(ID,1) - THETAW                                  30.82
          COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                    40.41
          AC2CEN = AC2(ID,IS,KCGRD(1))
!
          SWIND_EXP = 0.                                                  40.13
          SWIND_IMP = 0.                                                  40.13

          IF ( GROWW(ID,IS) ) THEN
!           *** term A ***
            IF ( SIGMA .GE. ( 0.7 * SIGPKD ) ) THEN
              ADUM = FACTA * (WIND10 * COSDIF)**4                         40.41
              ADUM = MAX ( 0. , ADUM / SIGTPI )
            ELSE
              ADUM = 0.
            END IF
!           *** term B; Note that BDUM is multiplied with a factor 5 ***
!           *** as in the DOLPHIN-B model                            ***
!
            BDUM = MAX( 0., ((WIND10 * CINV) * COSDIF-PWIND(3)))          40.41
            BDUM = FACTB * BDUM * 5.
            SWIND_EXP = ADUM + BDUM * AC2CEN                              40.13
!
          ELSE IF ( .NOT. GROWW(ID,IS) .AND. AC2CEN .GT. 0. ) THEN
!
!           *** for no energy dissipation outside the wind field     ***
!           *** TAUINV is set equal zero (as in the DOLPHIN-B model) ***
!
            IF ( COSDIF .LT. 0. ) THEN                                    40.41
              TAUINV = 0.
            ELSE
              TAUINV = ( SIGMA**2 * WIND10 * ABS(COSDIF) ) /              40.41
     &                 ( PWIND(4) * GRAV * TWOPI**2 )
            ENDIF
            SWIND_EXP = TAUINV * ALIMW(ID,IS)
            SWIND_IMP = -TAUINV
            ADUM = ALIMW(ID,IS)
            BDUM = TAUINV
          END IF
!
!         *** store results in IMATDA and IMATRA ***
!
          IMATRA(ID,IS) = IMATRA(ID,IS) + SWIND_EXP
          IMATDA(ID,IS) = IMATDA(ID,IS) - SWIND_IMP
          IF (TESTFL) PLWNDA(ID,IS,IPTST) = SWIND_EXP                     40.13
          IF (TESTFL) PLWNDB(ID,IS,IPTST) = SWIND_IMP                     40.13
          GENC0(ID,IS,1) = GENC0(ID,IS,1) + SWIND_EXP                     40.85
          GENC1(ID,IS,1) = GENC1(ID,IS,1) + SWIND_IMP                     40.85
!
!         *** test output ***
!
!         Value of ITEST changed from 10 to 110 to reduce test output     40.13
          IF ( TESTFL .AND. ITEST .GE. 110 ) THEN                         40.13
            WRITE(PRINTF,2004) IS, ID, GROWW(ID,IS), ADUM, BDUM           40.13
 2004       FORMAT(' WNDPAR: IS ID GROWW ADUM BDUM     :',                40.13
     &             2I3,2X,L1,2X,2E12.4)
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
!     Value of ITEST changed from 10 to 60 to reduce test output          40.13
      IF ( TESTFL .AND. ITEST .GE. 60 ) THEN                              40.13
        WRITE(PRINTF,*)
        WRITE(PRINTF,6051) IDWMIN, IDWMAX
 6051   FORMAT(' WNDPAR : IDWMIN IDWMAX     :',2I5)
        WRITE(PRINTF,6052) THETAW,WIND10,SIGPK,SIGPKD
 6052   FORMAT(' WNDPAR : Tw U10 Spk Spk,d   :',4E12.4)
        WRITE(PRINTF,7050) ETOTW, EDML, ALPM, ALPMD
 7050   FORMAT(' WNDPAR: ETOW EDML ALPM ALPMD:',4E12.4)
      ENDIF
!
      RETURN
!     end of subroutine WNDPAR
      END
!
!****************************************************************
!
      SUBROUTINE WINDP1 (WIND10     ,THETAW     ,
     &                   IDWMIN     ,IDWMAX     ,
     &                   FPM        ,UFRIC      ,
     &                   WX2        ,WY2        ,
     &                   ANYWND     ,SPCDIR     ,                         40.00
     &                   UX2        ,UY2        ,SPCSIG     )             30.70
!ADC     &                   UX2        ,UY2        ,SPCSIG     ,             41.20 30.70
!ADC     &                   NodeNumber )                                     41.20
!
!****************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
!ADC      USE WIND, ONLY: WindDrag                                            41.20
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
!     30.82: IJsbrand Haagsma
!     32.06: Roeland Ris
!     40.00: Nico Booij (Nonstationary boundary conditions)
!     40.41: Marcel Zijlema
!     41.20: Casey Dietrich
!
!  1. Updates
!
!     20.64,        : limited sector is now taken into account
!                     argument SPCDIR is added
!     30.70, Feb. 98: relative wind velocity is determined
!                     arguments UX2, UY2, SPCSIG added
!                     full common introduced
!     40.00, July 98: argument list changed: KCGRD removed
!     30.82, Oct. 98: Updated description of several variables
!     32.06, June 99: Reformulation of wind speed in terms of friction
!                     velocity for first and second generation
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.20, Jun. 10: use ADCIRC's new sector-based wind drag
!
!  2. Purpose
!
!        Computation of parameters derived from the wind for several
!        subroutines such as :
!                              SWIND1, SWIND2 SWIND3
!                              CUTOFF
!
!        Output of this subroutine :
!
!        WIND10 , THETAW, IDWMIN, IDWMAX , UFRIC, FPM
!
!     3. METHOD
!
!     a. For SWIND1 and SWIND2 :
!
!        SIGMA_FPM = 0.13 * GRAV * 2 * PI / WIND10
!
!     b. For SWIND3 (wind input according to Snyder (1981) ***
!
!       - wind friction velocity according to Wu (1982):
!          *                                                 -3
!         U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     c. For SWIND4 (wind input according to Janssen 1991)
!
!       - wind friction velocity:
!
!          UFRIC = sqrt ( CDRAG) * U10
!
!          for U10 < 7.5 m/s  ->  CDRAG = 1.2873.e-3
!
!          else wind friction velocity according to Wu (1982):
!
!          *                                                 -3
!         U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     d.
!        The Pierson Moskowitz radian frequency for a fully developed
!        sea state spectrum for all third generation wind input
!        models is equal to:
!
!                        grav
!        SIGMA_FPM  =  ---------
!                      28 UFRIC
!        -----------------------------------------------------------
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
!        IDWMIN           Minimum counter for spectral wind direction
!        IDWMAX           Maximum counter for spectral wind direction
!        IX               Counter of gridpoints in x-direction
!        IY               Counter of gridpoints in y-direction
!        MXC              Maximum counter of gridppoints in x-direction
!        MYC              Maximum counter of gridppoints in y-direction
!        KCGRD   int, i   Point index for grid point                      30.21
!        MCGRD   int, i   Maximum counter of gridpoints in space          30.21
!        ICMAX   int, i   Maximum counter for the points of the molecule  30.21
!
!        REALS:
!        ---------
!
!        PI          (3,14)
!        GRAV        Gravitational acceleration
!        THETAW      Mean direction of the relative wind vector
!        WIND10      Velocity of the relative wind vector
!        U10         wind velocity from SWANPREn
!        WDIC        mean wind direction from SWANPREn
!        FPM         PM frequency
!        UFRIC       Wind friction velocity
!
!        one and more dimensional arrays:
!        ---------------------------------
!        WX2,WY2     2D   Wind velocity array relative to a current
!        PWIND       1D   coefficients for wind expressions
!        ANYWND      1D   Indicator if wind input has to be taken into
!                         account for a bin
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
!
!     6. SUBROUTINES USED
!
!        NONE
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
!     ------------------------------------------------------------
!     If constant wind then
!       set parameters WIND10 and WDIC equal U10 and WDIC
!     else compute
!       wind velocity and mean wind direction
!     ----------------------------------------------------------
!     compute the minimum and maximum counters for wind source term
!     for which the wid input is active
!     compute the Pierson Moskowitz frequency for IWIND = 1, 2
!     compute the wind friction velocity and PM freq. for IWIND = 3
!     ------------------------------------------------------------
!     End of the subroutine WINDP1
!     ------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
!ADC      INTEGER, OPTIONAL :: NodeNumber                                     41.20
!ADC!
      INTEGER      IDWMIN ,IDWMAX                                         30.70
!
      REAL         WIND10 ,THETAW ,                                       30.70
     &             UFRIC  ,FPM    ,CDRAG  ,SDMEAN                         30.70
!
      REAL         WX2(MCGRD)   ,
     &             WY2(MCGRD)   ,
     &             UX2(MCGRD)   ,                                         30.70
     &             UY2(MCGRD)                                             30.70
!
      LOGICAL      ANYWND(MDC)
!
      REAL         AWX, AWY, RWX, RWY                                     30.70
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'WINDP1')
!
!     compute absolute wind velocity                                      30.70
      IF (VARWI) THEN
        AWX = WX2(KCGRD(1))
        AWY = WY2(KCGRD(1))
      ELSE
        AWX = U10 * COS(WDIC)
        AWY = U10 * SIN(WDIC)
      ENDIF
!     compute relative wind velocity                                      30.70
      IF (ICUR.EQ.0) THEN
        RWX = AWX
        RWY = AWY
      ELSE
        RWX = AWX - UX2(KCGRD(1))
        RWY = AWY - UY2(KCGRD(1))
      ENDIF
!     compute absolute value of relative wind velocity                    30.70
      WIND10 = SQRT(RWX**2+RWY**2)
      IF (WIND10.GT.0.) THEN
        THETAW = ATAN2 (RWY,RWX)
        THETAW = MOD ( (THETAW + PI2) , PI2 )
      ELSE
        THETAW = 0.
      ENDIF
!
!     *** compute the minimum and maximum counter for the active  ***
!     *** wind field :                                            ***
!     ***                                   .                     ***
!     *** IDWMAX =135 degrees             .<--mean wind direction ***
!     ***              \ o o o | o o o  .     (THETAW)            ***
!     ***                \ o o | o o  .  o                        ***
!     ***                  \ o | o  .  o o                        ***
!     ***                    \ |  .  o o o                        ***
!     ***      ----------------\--------------------              ***
!     ***                      | \ o o o o                        ***
!     ***                      |   \ o o o                        ***
!     ***                      |     \ o o                        ***
!     ***                             IDWMIN = 325 degrees        ***
!     ***                                                         ***
!
!     move ThetaW to the right interval, shifting + or - 2*PI             20.64
      SDMEAN = 0.5 * (SPCDIR(1,1) + SPCDIR(MDC,1))                        30.82
      IF (THETAW .LT. SDMEAN - PI) THETAW = THETAW + 2.*PI
      IF (THETAW .GT. SDMEAN + PI) THETAW = THETAW - 2.*PI
!
      IF ( (THETAW - 0.5 * PI) .LE. SPCDIR(1,1) ) THEN                    30.82
        IF ( (THETAW + 1.5 * PI) .GE. SPCDIR(MDC,1) ) THEN
          IDWMIN = 1
        ELSE
          IDWMIN = NINT ( (THETAW + 1.5*PI - SPCDIR(1,1)) / DDIR ) + 1    30.82
        ENDIF
      ELSE
        IDWMIN = NINT ( (THETAW - 0.5*PI - SPCDIR(1,1)) / DDIR ) + 1      30.82
      END IF
!
      IF ( (THETAW + 0.5 * PI) .GE. SPCDIR(MDC,1) ) THEN                  30.82
        IF ( (THETAW - 1.5 * PI) .LE. SPCDIR(1,1) ) THEN                  30.82
          IDWMAX = MDC
        ELSE
          IDWMAX = NINT ( (THETAW - 1.5 * PI - SPCDIR(1,1)) / DDIR ) + 1  30.82
        ENDIF
      ELSE
        IDWMAX = NINT ( (THETAW + 0.5 * PI - SPCDIR(1,1)) / DDIR ) + 1    30.82
      ENDIF
!
      IF ( IDWMIN .GT. IDWMAX) IDWMAX = MDC + IDWMAX
!
!     *** determine for which bin the wind input is active ***
!     *** initialize array for active wind input           ***
!
      DO ID = 1, MDC
        ANYWND(ID) = .FALSE.
      ENDDO
!
      IF ( TESTFL .AND. ITEST .GE. 30 ) THEN
          WRITE(PRINTF,500) IDWMIN, IDWMAX
 500      FORMAT(' WINDP1: IDWMIN IDWMAX :',2I15)
      ENDIF
!
      DO IDDUM = IDWMIN , IDWMAX
        ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
        ANYWND(ID) = .TRUE.
!
        IF ( TESTFL .AND. ITEST .GE. 40 ) THEN
          WRITE(PRINTF,400) IDDUM, ID, ANYWND(ID)
 400      FORMAT(' WINDP1: IDDUM ID ANYWND :',2I5,L4)
        ENDIF
!
      ENDDO
!
!     *** compute the Pierson Moskowitz frequency ***
!
      IF ( IWIND .EQ. 1 .OR. IWIND .EQ. 2 ) THEN
!
!       *** first and second generation wind wave model ***
!
        IF ( WIND10 .LT. PWIND(12) ) WIND10 = PWIND(12)
        FPM = 2. * PI * PWIND(13) * GRAV / WIND10
!
!       *** determine U friction in case predictor is obtained ***
!       *** with second genaration wave model                  ***
!
        IF ( WIND10 .GT. 7.5 ) THEN
          UFRIC = WIND10 * SQRT (( 0.8 + 0.065 * WIND10 ) * 0.001 )
        ELSE
          CDRAG = 0.0012873
          UFRIC = SQRT ( CDRAG ) * WIND10
        ENDIF
!
!     Reformulation of the wind speed in terms of friction velocity.      32.06
!     This formulation is based on Bouws (1986) and described in Delft    32.06
!     Hydraulics report H3515 (1999)                                      32.06
!
        WIND10 = WIND10 * SQRT(((0.8 + 0.065 * WIND10) * 0.001) /         32.06
     &                         ((0.8 + 0.065 * 15.   ) * 0.001))          32.06
!
      ELSE IF (IWIND .GE. 3 ) THEN
!
!       *** Calculate the wind friction velocity  ***
!       *** according to Wu (1982)                ***
!       *** apply cd-cap if appropriate           ***
!
        IF ( WIND10 .GT. 7.5 ) THEN
          CDRAG = ( 0.8 + 0.065 * WIND10 ) * 0.001
          CDRAG = MIN ( CDCAP, CDRAG )
!ADC          CDRAG = WindDrag(DBLE(WIND10),0.002D0,"Powell    ",NodeNumber)
        ELSE
          CDRAG = 0.0012873
        ENDIF
        UFRIC = SQRT ( CDRAG ) * WIND10
!
!       *** Wind friction velocity and PM-frequency ***
!
        UFRIC = MAX ( 1.E-15 , UFRIC)
        FPM =  GRAV / ( 28.0 * UFRIC )

      END IF
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 50 ) THEN
        WRITE(PRINTF,6050) KCGRD(1), MDC, MCGRD, IWIND                    30.21
 6050   FORMAT(' WINDP1:INDEX MDC MCGRD IWND:',4I5)
        WRITE(PRINTF,6052) THETAW,WIND10,WDIC,U10
 6052   FORMAT('       : THAW WIND10 WDIC U10 :',4E12.4)
        WRITE(PRINTF,6054) GRAV, PI, DDIR, VARWI
 6054   FORMAT('       : GRAV PI DDIR VARWI     :',3E12.4,L6)
        WRITE(PRINTF,6056) IDWMIN,IDWMAX,FPM, UFRIC
 6056   FORMAT('       : IDWMIN IDWMAX FPM UFR:',2I4,2E12.4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     end of subroutine WINDP1
      END
!
!****************************************************************
!
      SUBROUTINE WINDP2 (IDWMIN  ,IDWMAX  ,SIGPKD  ,FPM     ,
     &                   ETOTW   ,
     &                   AC2     ,SPCSIG  ,                               40.00
     &                   WIND10                                      )    30.70
!
!****************************************************************
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
!     30.72: IJsbrand Haagsma
!     30.70: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     20.72, Jan. 96: Integration modified, using FRINTF, FRINTH
!                     and PWTAIL(6)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: full common introduced, argument list changed
!                     ISFPM changed (in case of very high value of FPM)
!     40.00, July 98: argument list changed: KCGRD removed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Computation of wind sea energy spectrum for the second
!     generation wind growth model. Output of the subroutine
!     is ETOTW (total wind sea energy spectrum)
!
!  3. Method
!
!     Compute the wind sea energy spectrum : ETOTW
!
!             +90  inf
!             |   |
!     ETOTW = |   |        E(s,d) ds dd
!            -90  0.7*FPM
!
!     Compute the total energy density ETOTW for F > 0.7 FPM
!
!          ^  |
!       E()|  |          *
!             |        *   *
!             |              *
!             |       *      | *      / ETOTW(e)
!             |              | o  * /
!             |      *       | o o/o *
!             |              | o o o o o *
!             |     *        | o o o o o o o o*
!            0---------------|---------------------------
!             0          0.7*FPM              SIGMA --> s
!
!                   SIGMA MAX
!                  |
!      ETOTW(d) =  |  E(s,d)ds      ISFPM = FPM
!                 0.7 FPM
!
!      and over the interval +/- 90 degrees (according to
!      the mean wind direction as computed in WINDP1 )
!
!            IDWMAX = "135"
!               o
!      ETOTW =  o  ETOTW(d)dd
!
!         IDWMIN ="325"
!
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!        ISFPM       Counter in point just for the Pierson Moskowitz
!                    frequency
!        IDWMIN      Minimum counter for spectral wind direction
!        IDWMAX      Maximum counter for spectral wind direction
!        IX          Counter of gridpoints in x-direction
!        IY          Counter of gridpoints in y-direction
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MXC         Maximum counter of gridppoints in x-direction
!        MYC         Maximum counter of gridppoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!
!        DAC2        Difference in action density between two neighbour
!                    points (IS and IS+1)
!        DD          Directional band width
!        DSIG        Distance between FPM and next frequency
!        DUM_DS      Distance betwee the relative frequency and the
!                    next frequency value : LOGSIG(IS+1) - FPM
!        ETOT_       Dummy variables to compute the total energy
!        ETOTW       Total energy density over the sea spectrum, i.e.
!                    a > apm and a 180 degrees spectral interval
!        CNETOT      contribution of high frequency tail over the
!                    full spectrum
!        FPM         Pierson Moskowitz frequency
!        SIGFAC      Relative difference between two frequencies
!        THETAW      Mean direction of the relative wind vector
!        WIND10      Velocity of the relative wind vector
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU, SOURCE
!
!     6. SUBROUTINES USED
!
!        NONE
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
!     ------------------------------------------------------------
!     Determine the counter for the FPM frequency
!     integrate the wind sea energy spectrum
!     add the energy of the high frequency tail to the spectrum
!     ------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
!
      INTEGER  IDWMIN  ,IDWMAX  ,
     &         IDDUM   ,ID      ,IS      ,ISFPM
!
      REAL     ETOTW   ,FPM     ,SIG     ,ATOTD

      REAL     AC2(MDC,MSC,MCGRD)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'WINDP2')
!
!     *** compute wind sea energy spectrum for IS > 0.7 FPM       ***
!     *** minimum FPM is equal : 2 * pi * 0.13 * grav / pwind(12) ***
!     *** is equal 8 rad/s = 1.27 Hz                              ***
!
      ISFPM = MSC                                                         30.70
      FACINT = 0.                                                         30.70
      DO IS = 1, MSC
        SIG = SPCSIG(IS)                                                  30.72
        IF (FRINTH * SIG .GT. (0.7 * FPM) ) THEN
          ISFPM =  IS
          FACINT = (FRINTH - 0.7*FPM/SIG) / (FRINTH - 1./FRINTH)
          GOTO 11
        END IF
      ENDDO
 11   CONTINUE
!
!     *** calculate the energy in the wind sea part of the spectrum ***
!     *** from ISFPM.                                               ***
!
      ETOTW = 0.
      DO IS = ISFPM, MSC
        SIG = SPCSIG(IS)                                                  30.72
        ATOTD = 0.
        DO IDDUM = IDWMIN, IDWMAX
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          ATOTD = ATOTD + AC2(ID,IS,KCGRD(1))                             30.21
        ENDDO
        IF (IS.EQ.ISFPM) THEN
          ETOTW = ETOTW + FACINT * FRINTF * SIG**2 * DDIR * ATOTD
        ELSE
          ETOTW = ETOTW + FRINTF * SIG**2 * DDIR * ATOTD
        ENDIF
      ENDDO
!     add high-frequency tail:
      ETOTW = ETOTW + PWTAIL(6) * SIG**2 * DDIR * ATOTD
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 70 ) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,6050) IWIND,IDWMIN,IDWMAX, ISFPM, ETOTW
 6050   FORMAT(' WINDP2: IWND IDWMIN IDWMAX ISFPM ETOTW:',4I6,1X,E12.4)
      END IF
!
      RETURN
!     end of subroutine WINDP2
      END
!
!********************************************************************
!
      SUBROUTINE WINDP3 (ISSTOP  ,ALIMW   ,AC2     ,
     &                   GROWW   ,IDCMIN  ,IDCMAX  )                      40.41
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     1. UPDATE
!
!        40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!     2. PURPOSE
!
!        Reduce the energy density in spectral direction direct after
!        solving the tri-diagonal matrix if the energy density level is
!        larger then the upper bound limit given by a Pierson Moskowitz
!        spectrum. This is only carried out if a particular wave
!        component is 'growing'.
!
!        If the energy density in a bin is larger than the upper bound
!        limit (for instance when crossing wind seas are present) then
!        the energy density level is a lower limit.
!
!     3. METHOD
!
!        The upper bound limit is given by:
!
!                       2                  -4
!               ALPHA  g              sigma     2     2
!    A(s,t) = ----------- exp ( -5/4 (----)   ) --- cos ( d - dw )
!               sigma^6                FPM       pi
!
!         in which ALPHA is wind sea dependent (see subroutine
!         SWIND2) :
!
!                      /                                 C60   \
!                     |                / E_windsea g^2 \        |
!         ALPHA = MIN | 0.0081 ,  C50 |  ------------   |       |
!                      \               \    U_10^4     /        /
!
!
!     4. PARAMETERLIST
!
!        INTEGERS :
!        ----------
!        IX IY            Grid point in geographical space
!        MDC ,MSC         Counters in spectral space
!        ISSTOP           Maximum frequency that fall within a sweep
!
!        ARRAYS:
!        -------
!        AC2      4D      Action density as funciton of x,y,s,t
!        ALIMW    2D      Contains the action density upper bound
!                         limit regardind spectral action density per
!                         spectral bin (A(s,t))
!        GROWW    2D      Logical array which determines is there is
!                         a) generation  ( E < E_lim -> .TRUE.  ) or
!                         b) dissipation ( E > E_lim -> .FALSE. )
!        IDCMIN   1D      Frequency dependent minimum counter
!        IDCMAX   1D      Frequency dependent minimum counter
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
!
!     6. SUBROUTINES USED
!
!        NONE
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every spectral bin
!     If wave input for a bin is true (GROWW(ID,IS) = .TRUE.) and
!        if wave action is larger then maximum then reduce
!        wave action to limit its value to upper boundary limit
!     else if no growth (see subroutine SWIND1 and SWIND2), check
!       if the energy density level in the wind sea part
!       of the spectrum is lower than upper bound limit. If so
!       set energy level equal to upper bound limit
!     endif
!   ------------------------------------------------------------
!   End of the subroutine WINDP3
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER     IS      ,ID      ,ISSTOP  ,IDDUM
!
      INTEGER     IDCMIN(MSC)       ,
     &            IDCMAX(MSC)
!
      REAL        AC2CEN
!
      REAL        AC2(MDC,MSC,MCGRD)    ,
     &            ALIMW(MDC,MSC)
!
      LOGICAL     GROWW(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'WINDP3')
!
!     *** limit the action density spectrum ***
!
      DO IS = 1, ISSTOP
        DO IDDUM = IDCMIN(IS) , IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          AC2CEN = AC2(ID,IS,KCGRD(1))
          IF ( GROWW(ID,IS) .AND. AC2CEN .GT. ALIMW(ID,IS) )
     &      AC2(ID,IS,KCGRD(1)) = ALIMW(ID,IS)
          IF ( .NOT. GROWW(ID,IS) .AND. AC2CEN .LT. ALIMW(ID,IS) )
     &      AC2(ID,IS,KCGRD(1)) = ALIMW(ID,IS)
!
          IF (TESTFL .AND. ITEST .GE. 50) THEN
             WRITE(PRINTF,300) IS,ID,GROWW(ID,IS),AC2CEN,ALIMW(ID,IS)
 300         FORMAT(' WINDP3 : IS ID GROWW AC2CEN ALIM:',2I4,L4,2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (TESTFL .AND. ITEST .GE. 50) THEN
        WRITE(PRINTF,4000) KCGRD(1),ISSTOP,MSC,MDC,MCGRD
 4000   FORMAT(' WINDP3 : POINT ISSTOP MSC MDC MCGRD :',5I5)
      END IF
!
      RETURN
!     end of subroutine WINDP3
      END
!
!****************************************************************
!
      SUBROUTINE SWIND0 (IDCMIN  ,IDCMAX  ,ISSTOP  ,
     &                   SPCSIG  ,THETAW  ,ANYWND  ,
     &                   UFRIC   ,FPM     ,PLWNDA  ,
     &                   IMATRA  ,SPCDIR  ,GENC0   )
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.41: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     40.41, Aug. 04: COS(THETA-THETAW) replaced by sumrule to make it cheaper
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.85, Aug. 08: store wind input for output purposes
!
!  2. PURPOSE
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Linear wind input term according to Cavaleri and
!         Malanotte-Rizzoli (1981)
!
!  3. METHOD
!
!     To ensure wave growth when no wave energy is present in the
!     numerical model a linear growth term is used (see
!     Cavaleri and Malanotte-Rizzoli 1981). Contributions for
!     frequencies lower then FPM have been eliminated buy aa filter :
!
!                  -3
!            1.5*10      *                       4      sigma -4
!  A = Sin = ------- { U  max[0, (cos(d - dw )] } exp{-(-----)  } / Jac
!               2                                        FPM
!              g
!
!        With Jac = Jacobian =  2 pi sigma
!
!        The Pierson Moskowitz radian frequency for a fully developed
!        sea state spectrum is as follows (computed in WINDP2 :
!
!                1       g
!        FPM  = ---- ---------  * 2 pi
!               2 pi  28 UFRIC
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
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!        IDWMIN      Minimum counter for spectral wind direction
!        IDWMAX      Maximum counter for spectral wind direction
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!        FPM         Pierson Moskowitz frequecy (WAM)
!        GRAV        Gravitatiuonal acceleration
!        SWINEO      Coefficient stored in IMATRA
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        UFRIC       Wind friction velocity
!
!        one and more dimensional arrays:
!        ---------------------------------
!        IMATRA    2D    Coefficients of right hand side of matrix
!        ANYWND    1D    Determines the number of directional bins
!                        that fall within the wind sea part of the
!                        spectrum ( thetaw +- 90 degrees)
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
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
!        ---
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every spectral bin that faal witin the sweep considered:
!     compute the source term Sinl and store the results in array
!   --------------------------------------------------------------
!   End of the subroutine SWIND0
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER  IDDUM   ,ID      ,IS      ,ISSTOP
!
      REAL     FPM     ,UFRIC   ,THETA   ,THETAW  ,
     &         SWINEA  ,SIGMA   ,TEMP1   ,TEMP2   ,
     &         CTW     ,STW     ,COSDIF  ,                                40.41
     &         TEMP3   ,FILTER
!
      REAL    IMATRA(MDC,MSC)      ,
     &        PLWNDA(MDC,MSC,NPTST)                                       40.00
      REAL    GENC0(MDC,MSC,MGENR)                                        40.85
!
      INTEGER IDCMIN(MSC)          ,
     &        IDCMAX(MSC)
!
      LOGICAL ANYWND(MDC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWIND0')
!
!     *** calculate linear wind input term ***
!
      CTW = COS(THETAW)                                                   40.41
      STW = SIN(THETAW)                                                   40.41
      FPM =  GRAV / ( 28.0 * UFRIC )
      TEMP1 = PWIND(31) / ( GRAV**2 * 2. * PI )                           7/MAR
      DO IS = 1, ISSTOP
        SIGMA  = SPCSIG(IS)                                               30.72
!
!       ****            ARGU   =  FPM / SIGMA                     ***
!       **** the value of ARGU was change for MIN () because for  ***
!       **** values of fpm/sigma too small could be some problems ***
!       **** with some computers to handle small numbers          ***
        ARGU   = MIN (2., FPM / SIGMA)                                    30.00
        FILTER = EXP ( - ARGU**4 )                                        20.87
        TEMP2  = TEMP1 / SIGMA
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          IF ( ANYWND(ID) .AND. SIGMA .GE. (0.7 * FPM) ) THEN
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            TEMP3  = ( UFRIC *  MAX( 0. , COSDIF))**4                     40.41
            SWINEA = MAX( 0. , TEMP2 * TEMP3 * FILTER )
            IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEA
            IF(TESTFL) PLWNDA(ID,IS,IPTST) = SWINEA/SIGMA                 40.85 40.00
            GENC0(ID,IS,1) = GENC0(ID,IS,1) + SWINEA/SIGMA                40.85
!
!           *** test output ***
!
            IF (ITEST .GE. 80 .AND. TESTFL )
     &      WRITE (PRTEST, 333) ID, IS, FILTER, SWINEA
 333        FORMAT (' ID IS FILTER  WIND SOURCE (IMATRA)',
     &               2I4, 1X, 2(1X,E11.4))
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST.GE. 60.AND.TESTFL) THEN
        WRITE(PRINTF,400) KCGRD(1), THETAW*180./PI
 400    FORMAT(' SWIND0: POINT  THETAW       :',I5,E12.4)
        WRITE(PRINTF,500) TEMP1, FPM, UFRIC
 500    FORMAT(' SWIND0: TEMP1 FPM UFRC     :',3E12.4)
        WRITE(PRINTF,*)
        IF (ITEST.GE. 120.AND.TESTFL) THEN                                 24/MAR
          DO IS = 1, ISSTOP
            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
              ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
              WRITE(PRINTF,100) IS,ID,ANYWND(ID)
 100          FORMAT(' IS ID ANYWND : ', 2I5,1X,L2)
            ENDDO
          ENDDO
        ENDIF
      END IF
!
      RETURN
!     end of subroutine SWIND0
      END
!
!****************************************************************
!
      SUBROUTINE SWIND3 (SPCSIG  ,THETAW  ,
     &                   KWAVE   ,IMATRA  ,GENC0   ,
     &                   IDCMIN  ,IDCMAX  ,AC2     ,UFRIC   ,
     &                   FPM     ,PLWNDA  ,ISSTOP  ,SPCDIR  ,ANYWND)
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.41: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     40.41, Aug. 04: COS(THETA-THETAW) replaced by sumrule to make it cheaper
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.85, Aug. 08: store wind input for output purposes
!
!  2. PURPOSE
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Exponential input term, (Snyder et al. 1981, which
!         expression has been modified by Komen et al. 1984).
!
!         This input term should be combinated with the dissipation
!         term of Komen et al. (1984)
!
!  3. METHOD
!
!     The exponential term used is taken from Snyder et al. (1981)
!     and Komen et al. (1984):
!
!     Sin (s,d) =  B*E(s,d)
!        e                             *
!                                 28 U cos( d - dw )
!     B = max(0. , (0.25 rhoaw( -------------------  -1 ) )) sigma
!                                   sigma / kwave
!
!     with :
!
!        *                                                -3
!      U  =  UFRIC = wind10 sqrt( (0.8 + 0.065 wind10 ) 10  )
!
!     UFRIC is computed in WINDP1
!
!
!     The Pierson Moskowitz radian frequency for a fully developed
!     sea state spectrum is as follows (computed in WINDP1):
!
!             1       g
!     FPM  = ---- ---------  * 2 pi
!            2 pi  28 UFRIC
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
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!        FPM         Pierson Moskowitz frequecy (WAM)
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        UFRIC       Wind friction velocity
!
!        one and more dimensional arrays:
!        ---------------------------------
!        KWAVE     2D    Wavenumber
!        LOGSIG    1D    Logaritmic distribution of frequency
!        IMATRA    2D    Coefficients of right hand side of matrix
!        PWIND     1D    Wind coefficients
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!        ANYWND    1D    Wind input for bin considered
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
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
!        ---
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every spectral bin that fall in teh sweep considered
!     compute the source term and store the results in IMATRA
!   --------------------------------------------------------------
!   End of the subroutine SWIND3
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER  IDDUM ,ID    ,IS    ,ISSTOP
!
      REAL     FPM   ,UFRIC ,THETA ,THETAW,SIGMA ,SWINEB,TEMP1,
     &         CTW   ,STW   ,COSDIF,                                      40.41
     &         TEMP2 ,TEMP3 ,CINV
!
      REAL    AC2(MDC,MSC,MCGRD)   ,
     &        IMATRA(MDC,MSC)      ,
     &        KWAVE(MSC,ICMAX)     ,
     &        PLWNDA(MDC,MSC,NPTST)                                       40.00
      REAL  :: GENC0(MDC,MSC,MGENR)                                       40.85
!
      INTEGER IDCMIN(MSC)          ,
     &        IDCMAX(MSC)
!
      LOGICAL  ANYWND(MDC)
!
!/T      LOGICAL  IMP_EXP
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWIND3')
!
      CTW   = COS(THETAW)                                                 40.41
      STW   = SIN(THETAW)                                                 40.41
      TEMP1 = 0.25 * PWIND(9)
      TEMP2 = 28.0 * UFRIC
      DO IS = 1, ISSTOP
        SIGMA = SPCSIG(IS)                                                30.72
        CINV  = KWAVE(IS,1) / SIGMA
        TEMP3 = TEMP2 * CINV
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          IF ( ANYWND(ID) ) THEN
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            SWINEB = TEMP1 * ( TEMP3 * COSDIF - 1.0 )                     40.41
            SWINEB = MAX ( 0. , SWINEB * SIGMA )
!
            IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,KCGRD(1))
            IF (TESTFL) PLWNDA(ID,IS,IPTST) = SWINEB*AC2(ID,IS,KCGRD(1))  40.85 40.00
            GENC0(ID,IS,1) = GENC0(ID,IS,1) + SWINEB*AC2(ID,IS,KCGRD(1))  40.85
!
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST.GE. 80.AND.TESTFL) THEN                                   40.00
        WRITE(PRTEST,6000) KCGRD(1), THETAW*180./PI
 6000   FORMAT(' SWIND3: POINT  THETAW        :',I5,E12.4)
        WRITE(PRTEST,6100) TEMP1, FPM, UFRIC
 6100   FORMAT(' SWIND3: TEMP1 FPM UFRC     :',3E12.4, /,
     &         '  IS ID1 ID2       Wind source term')
        DO IS = 1, MSC
          WRITE(PRTEST,6200) IS, IDCMIN(IS), IDCMAX(IS),
     &    (PLWNDA(ID,IS,IPTST), ID=IDCMIN(IS), IDCMAX(IS))
 6200     FORMAT(3I4, 600e12.4)
        ENDDO
        WRITE(PRTEST,*)
      END IF
!
      RETURN
!     end of subroutine SWIND3
      END
!
!****************************************************************
!
      SUBROUTINE SWIND4 (IDWMIN  ,IDWMAX  ,
     &                   SPCSIG  ,WIND10  ,THETAW  ,XIS     ,
     &                   DD      ,KWAVE   ,IMATRA  ,GENC0   ,
     &                   IDCMIN  ,IDCMAX  ,AC2     ,UFRIC   ,
     &                   PLWNDA  ,ISSTOP  ,ITER    ,USTAR   ,ZELEN   ,
     &                   SPCDIR  ,ANYWND  ,IT               )
!
!******************************************************************
!
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!     40.61: Roop Lalbeharry
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     40.02, Oct. 00: References to CDRAGP and TAUWP removed
!     40.31, Jul. 03: correction calculation TAUDIR in test output
!     40.41, Aug. 04: COS(THETA-THETAW) replaced by sumrule to make it cheaper
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Nov. 06: improvements to WAM4 based on WAM4.5
!     40.85, Aug. 08: store wind input for output purposes
!
!  2. Purpose
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     1)  Computation of the exponential input term based on a
!         quasi linear theory developped by Janssen (1989, 1991).
!         This formulation should be used in combination with the
!         whitecapping dissipation source term according to
!         Janssen (1991) and Mastenbroek et al. (1993)
!
!  3. Method
!
!     The exponential term for the wind input used is taken
!     from Janssen (1991):
!
!     Sin (s,d) =  B*E(s,d)
!        e
!                                        *
!                               / kwave U  \    2
!     B = max(0. , (beta rhoaw | --------- | cos( d - dw ) sigma))
!                               \ sigma   /
!
!     The friction velocity is a fucntion of the roughness
!     length Ze. A first gues for U* is given by Wu (1982),
!     which is computed in subroutine WINDP2.
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
!        IDWMIN      Minimum counter for spectral wind direction
!        IDWMAX      Maximum counter for spectral wind direction
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!        DD          Directional band width
!        GRAV        Gravitational acceleration
!        SWINEO      Coefficient stored in IMATRA
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        WIND10      Velocity of the relative wind vector
!        UFRIC       Wind friction velocity
!        ZALP        Wave growth parameter used in WAM4                   40.61
!
!        one and more dimensional arrays:
!        ---------------------------------
!        KWAVE     2D    Wavenumber
!        IMATRA    2D    Coefficients of right hand side of matrix
!        PWIND     1D    Wind coefficients
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!        USTAR     2D    Friction velocity at previous iteration level
!        ZELEN     2D    Roughness length at previous iteration level
!        ANYWND    1D    Determine if wind input is active for bin
!
!        PWIND(9)  1D    Rho air / rho water
!        ------------
!        PWIND(14) 1D    alfa (which is tuned at 0.01)
!        PWIND(15) 1D    Kappa ( 0.41)
!        PWIND(16) 1D    Rho air
!        PWIND(17) 1D    Rho water
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
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
!        ---
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For
!   --------------------------------------------------------------
!   End of the subroutine SWIND4
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!
!***********************************************************************
!
      INTEGER  IDWMAX  ,IDWMIN  ,IDDUM   ,ID      ,ISSTOP  ,IS
!
      REAL     THETA  ,THETAW ,DD     ,SWINEB ,WIND10 ,
     &         ZO     ,ZE     ,BETA1  ,BETA2  ,UFRIC  ,UFRIC2 ,DS     ,
     &         ZARG   ,ZLOG1  ,ZLOG2  ,ZCN1   ,ZCN2   ,ZCN    ,XIS    ,
     &         SIGMA  ,SIGMA1 ,SIGMA2 ,WAVEN  ,WAVEN1 ,WAVEN2 ,TAUW   ,
     &         TAUTOT ,TAUDIR ,COS1   ,COS2   ,CW1    ,RHOA   ,RHOW   ,
     &         RHOAW  ,ALPHA  ,XKAPPA ,F1     ,TAUWX  ,TAUWY  ,SE1    ,
     &         CTW    ,STW    ,COSDIF ,                                   40.41
     &         SE2    ,SINWAV ,COSWAV
      REAL     ZALP                                                       40.61
!
      REAL    AC2(MDC,MSC,MCGRD)   ,                                      30.21
     &        IMATRA(MDC,MSC)      ,
     &        KWAVE(MSC,ICMAX)     ,
     &        PLWNDA(MDC,MSC,NPTST),
     &        USTAR(MCGRD)         ,
     &        ZELEN(MCGRD)
      REAL    GENC0(MDC,MSC,MGENR)                                        40.85
!
      INTEGER IDCMIN(MSC)          ,
     &        IDCMAX(MSC)
!
      LOGICAL ANYWND(MDC)                                                 40.41
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWIND4')
!
!
!     *** initialization ***
!
      ALPHA  = PWIND(14)
      XKAPPA = PWIND(15)
      RHOA   = PWIND(16)
      RHOW   = PWIND(17)
      RHOAW  = RHOA / RHOW
      ZTEN   = 10.
      RATIO  = 0.75
      BETAMX = 1.2
      F1     = BETAMX / XKAPPA**2
      ZALP   = 0.011                                                      40.61
      CTW    = COS(THETAW)                                                40.41
      STW    = SIN(THETAW)                                                40.41
!
      IF ( NSTATC.EQ.1 .AND. IT.EQ.1 ) THEN                               40.41 40.00
!
!        *** nonstationary and first time step (the number of        ***
!        *** iterations however still can increase per time step     ***
!
        ZO     = ALPHA * UFRIC * UFRIC / GRAV
        ZE     = ZO / SQRT( 1. - RATIO )
        USTAR(KCGRD(1)) = UFRIC                                           30.21
        ZELEN(KCGRD(1)) = ZE                                              30.21
      ELSE IF ( NSTATC.EQ.0 .AND. ICOND.EQ.4 .AND. ITER .EQ. 1 ) THEN     40.41 40.00
!
!        *** non-first stationary computations and first iteration   ***
!
        ZO     = ALPHA * UFRIC * UFRIC / GRAV
        ZE     = ZO / SQRT( 1. - RATIO )
        USTAR(KCGRD(1)) = UFRIC                                           30.21
        ZELEN(KCGRD(1)) = ZE                                              30.21
      ELSE IF ( NSTATC.EQ.0 .AND. ICOND.NE.4 .AND. ITER .EQ. 2 ) THEN     40.41 40.00
!
!        *** first stationary computation (this subroutine is never ***
!        *** excecuted anyway, this subroutine in entered after 1   ***
!        *** iteration) and thus calculate ZO and ZE as a first     ***
!        *** prediction only and only in the second sweep           ***
!
        ZO     = ALPHA * UFRIC * UFRIC / GRAV
        ZE     = ZO / SQRT( 1. - RATIO )
        USTAR(KCGRD(1)) = UFRIC
        ZELEN(KCGRD(1)) = ZE
      ELSE
!
!       *** calculate wave stress using the value of the  ***
!       *** velocity U* and roughness length Ze from the  ***
!       *** previous iteration                            ***
!
        UFRIC = USTAR(KCGRD(1))
        ZE    = ZELEN(KCGRD(1))
!
        TAUW   = 0.
        TAUWX  = 0.
        TAUWY  = 0.
        TXHFR  = 0.
        TYHFR  = 0.
!
!       *** use old friction velocity to calculate wave stress ***
!
        UFRIC2 = UFRIC * UFRIC
!
        DO IS = 1, MSC-1
          SIGMA1 = SPCSIG(IS)                                             30.72
          SIGMA2 = SPCSIG(IS+1)                                           30.72
          WAVEN1 = KWAVE(IS,1)
          WAVEN2 = KWAVE(IS+1,1)
          DS     = SIGMA2 - SIGMA1
          CW1    = SIGMA1 / WAVEN1
          CW2    = SIGMA2 / WAVEN2
          ZCN1   = ALOG ( GRAV * ZE / CW1**2 )
          ZCN2   = ALOG ( GRAV * ZE / CW2**2 )
          X1     = (UFRIC/CW1 + ZALP)**2                                  40.61
          X2     = (UFRIC/CW2 + ZALP)**2                                  40.61
          DO IDDUM = IDWMIN, IDWMAX
            ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            SINWAV = SPCDIR(ID,3)                                         40.41
            COSWAV = SPCDIR(ID,2)                                         40.41
            COS1   = MAX ( 0. , COSDIF )                                  40.41
            COS2   = COS1 * COS1
            BETA1  = 0.
            BETA2  = 0.
!
!           *** Miles constant Beta ***
!
            IF ( COS1 .GT. 0.01 ) THEN
              ZARG1 = XKAPPA / ( (UFRIC / CW1 + ZALP ) * COS1 )           40.61
              ZARG2 = XKAPPA / ( (UFRIC / CW2 + ZALP ) * COS1 )           40.61
              ZLOG1 = ZCN1 + ZARG1
              ZLOG2 = ZCN2 + ZARG2
              IF (ZLOG1.LT.0.) BETA1 = F1 * X1 * EXP (ZLOG1) * ZLOG1**4   40.61
              IF (ZLOG2.LT.0.) BETA2 = F1 * X2 * EXP (ZLOG2) * ZLOG2**4   40.61
            ENDIF
!
!           *** calculate wave stress by integrating input source ***
!           *** term in x- and y direction respectively           ***
!
            SE1 = BETA1 * SIGMA1**3 * AC2(ID,IS  ,KCGRD(1))               40.61
            SE2 = BETA2 * SIGMA2**3 * AC2(ID,IS+1,KCGRD(1))               40.61
!
            TAUWX = TAUWX + 0.5 * ( SE1 + SE2 ) * DS * COSWAV * COS2
            TAUWY = TAUWY + 0.5 * ( SE1 + SE2 ) * DS * SINWAV * COS2
!
!           *** test output ***
!
            IF (ITEST.GE. 40 .AND. TESTFL) THEN
              WRITE(PRINTF,105) IS, ID, UFRIC, ZE
  105         FORMAT(' SW4: IS ID UFRIC ZE     :',2I4,2E12.4)
              WRITE(PRINTF,106) ZLOG1, ZLOG2, BETA1, BETA2
  106         FORMAT(' SW4: ZOLG1-2 BETA1 BETA2:',4E12.4)
              IF (ABS(TAUWX).GT.0. .OR. ABS(TAUWY).GT.0.) THEN            40.31
                TAUDIR = ATAN2 ( TAUWX, TAUWY )                           40.31
              ELSE                                                        40.31
                TAUDIR = 0.                                               40.31
              ENDIF                                                       40.31
              TAUDIR = MOD ( (TAUDIR + 2. * PI) , (2. * PI) )             40.31
              WRITE(PRINTF,107) TAUWX, TAUWY, TAUDIR*180./PI
  107         FORMAT(' SW4: TAUWX TAUWY TAUDIR :',3E12.4)
            ENDIF
!
          ENDDO
        ENDDO
!
!       *** determine effect of high frequency tail to wave stress ***
!       *** assuming deep water conditions                         ***
!
        GAMHF =  XKAPPA * GRAV / UFRIC
        SIGMAX = SPCSIG(MSC)                                              30.72
        SIGHF1 = SIGMAX
        DO J=1, 50
          SIGHF2 = XIS * SIGHF1
          DS     = SIGHF2 - SIGHF1
          ZCNHF1 = ALOG ( ZE * SIGHF1**2 / GRAV )
          ZCNHF2 = ALOG ( ZE * SIGHF2**2 / GRAV )
          AUX    = 0.0
          SIGHF1 = SIGMAX                                                 40.61
          CW1    = GRAV/SIGHF1                                            40.61
          CW2    = GRAV/SIGHF2                                            40.61
          DO IDDUM = IDWMIN, IDWMAX
            ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            SINWAV = SPCDIR(ID,3)                                         40.41
            COSWAV = SPCDIR(ID,2)                                         40.41
            COS1   = MAX ( 0. , COSDIF )                                  40.41
            COS2   = COS1 * COS1
            COS3   = COS2 * COS1                                          40.61
            BETA1  = 0.0
            BETA2  = 0.0
!
            IF ( COS1 .GT. 0.01 ) THEN
!             *** beta is independent of direction ! ***
              ZAHF1 = XKAPPA / ( UFRIC / CW1 + ZALP )                     40.61
              ZAHF2 = XKAPPA / ( UFRIC / CW2 + ZALP )                     40.61
              ZLOG1 = ZCNHF1 + ZAHF1
              ZLOG2 = ZCNHF2 + ZAHF2
              IF ( ZLOG1 .LT. 0. ) BETA1 = F1 * EXP (ZLOG1) * ZLOG1**4
              IF ( ZLOG2 .LT. 0. ) BETA2 = F1 * EXP (ZLOG2) * ZLOG2**4
              AUX = AUX + BETA1 + BETA2
            ENDIF
!
!           *** calculate contribution of high frequency tail to ***
!           *** wave stress by integrating input source term in  ***
!           *** x- and y direction respectively                  ***
!
            FACHFR = SIGMAX**6 * AC2(ID,MSC,KCGRD(1)) * COS2 / GRAV**2    30.21

            SE1 = FACHFR * BETA1 / SIGHF1
            SE2 = FACHFR * BETA2 / SIGHF2
!
            TXHFR = TXHFR + 0.5 * ( SE1 + SE2 ) * DS * COSWAV
            TYHFR = TYHFR + 0.5 * ( SE1 + SE2 ) * DS * SINWAV
!
!           *** if coeffcient BETA = 0. for a frequency over ***
!           *** all directions is zero skip loop             ***
!
            IF ( AUX .EQ. 0. ) GOTO 5000
!
          ENDDO
!
          IF (ITEST.GE. 45 ) THEN
            WRITE(PRINTF,407) XIS, SIGHF1, SIGHF2, J
  407       FORMAT(' SW4: XIS SIGHF1 SIGHF2 J :',3E12.4,I4)
            WRITE(PRINTF,437) TXHFR, TYHFR, BETA1, BETA2
  437       FORMAT(' SW4: TXHFR TYHFR BETA1,2 :',4E12.4)
          ENDIF
!
          SIGHF1 = SIGHF2
        ENDDO
 5000   CONTINUE
!
        IF ( ITEST .GE. 45 ) THEN
          WRITE(PRINTF,321) TAUWX, TAUWY, TXHFR, TYHFR
 321      FORMAT(' SW4: Twx Twy Thfx Thfy:',4E12.4)
        ENDIF
!
        TAUTOT = RHOA * UFRIC2
!       *** wave stress ***
        TAUWX  = TAUWX + TXHFR * UFRIC2                                   40.61
        TAUWY  = TAUWY + TYHFR * UFRIC2                                   40.61
        IF (ABS(TAUWX).GT.0. .OR. ABS(TAUWY).GT.0.) THEN
          TAUDIR = ATAN2 ( TAUWX, TAUWY )
        ELSE
          TAUDIR = 0.
        ENDIF
        TAUDIR = MOD ( (TAUDIR + 2. * PI) , (2. * PI) )
        TAUW   = RHOA * DD * SQRT ( TAUWX**2 + TAUWY**2 )                 40.61
        TAUW   = MIN ( TAUW , 0.999 * TAUTOT )
!
        IF ( ITEST .GE. 45 ) THEN
          RATIO = TAUW / TAUTOT
          WRITE(PRINTF,301) TAUW, TAUTOT, RATIO, KCGRD(1)                 30.21
 301      FORMAT(' SW4: Tauw Taut  ratio :',3E12.4,' in ',I5)
        ENDIF
!
        DO II = 1, 20
!         *** start iteration process ***
          FA = SQRT ( 1. - TAUW / TAUTOT )
          FB = ZTEN * RHOA * GRAV / ALPHA
          FC = FA * ( FB / TAUTOT  - 1. )
          FD = SQRT ( TAUTOT )
          FE = ALOG ( FC + 1. )
!
!         *** calculate function value and derivative in ***
!         *** numerical point considered                 ***
!
          FCEN = FD * FE - SQRT(RHOA) * WIND10 * XKAPPA
          FF1  = 0.5 * FE / FD
          FF2  = 0.5 * TAUW * FC / FA**2 - FA * FB
          FF3  = TAUTOT**1.5 * ( FC + 1. )
          DCEN = FF1 + FF2 / FF3
!
!         *** new total stress ***
!
          TAUNEW = TAUTOT - FCEN / DCEN
!
          IF ( ITEST .GT. 30 .AND. TESTFL ) THEN
            WRITE(PRINTF,440) TAUTOT, TAUNEW, FCEN, DCEN, II
 440        FORMAT(' SW4: Tt Tnew Fcn DFcn II:',4E12.4,I2)
            WRITE(PRINTF,450) FA, FB, FC, FD
 450        FORMAT(' SW4: FA FB FC FD        :',4E12.4)
            WRITE(PRINTF,460) FE, FF1, FF2, FF3
 460        FORMAT(' SW4: FE FF1 FF2 FF3     :',4E12.4)
          ENDIF
!
          IF ( TAUNEW .LE. TAUW ) TAUNEW = .5 * (TAUTOT + TAUW)           20.81
          IF ( ABS ( TAUNEW - TAUTOT ) .LE. 1.E-5 ) GOTO 3000
!
          TAUTOT = TAUNEW
        ENDDO
 3000   CONTINUE
!
        UFRIC  = SQRT ( TAUTOT / RHOA )
!
        IF ( ITEST .GE. 20 .AND. TESTFL ) THEN
          WRITE(PRINTF,200) KCGRD(1)
 200      FORMAT(' SW4: Values after Newton-Raphson in point:',I5)
          WRITE(PRINTF,206) TAUW, TAUTOT, TAUW/TAUTOT, UFRIC
 206      FORMAT(' SW4: Tauw Taut rat Us :',4E12.4)
          WRITE(PRINTF,*)
        ENDIF
!
        ZO     = ALPHA * UFRIC * UFRIC / GRAV
        ZE     = ZO / SQRT ( 1. - TAUW / TAUTOT )
!
        USTAR(KCGRD(1)) = UFRIC
        ZELEN(KCGRD(1)) = ZE
!
      ENDIF
!
! ----->
!
!     *** calculate critical height and Miles parameter and  ***
!     *** calculate input source term B for with the updated ***
!     *** values of UFRIC and ZE                             ***
!
      UFRIC2 = UFRIC * UFRIC
!
      DO IS = 1, ISSTOP
        SIGMA  = SPCSIG(IS)                                               30.72
        WAVEN  = KWAVE(IS,1)
        CW1    = SIGMA / WAVEN
        ZCN    = ALOG ( GRAV * ZE / CW1**2 )
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          IF ( ANYWND(ID) )  THEN
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            COS1   = MAX ( 0. , COSDIF )                                  40.41
            COS2   = COS1 * COS1
            XFAC2  = ( UFRIC / CW1 + ZALP )**2                            40.61
            BETA   = 0.
            IF ( COS1 .GT. 0.01 ) THEN
              ZARG = XKAPPA / ( (UFRIC / CW1 + ZALP ) * COS1 )            40.61
              ZLOG = ZCN + ZARG
              IF ( ZLOG .LT. 0. ) BETA = F1 * EXP (ZLOG) * ZLOG**4
            ENDIF
!
!           *** compute the factor B and store result in array ***
!
            SWINEB = RHOAW * BETA * XFAC2 * COS2 * SIGMA
            IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,KCGRD(1))  30.21
            IF (TESTFL) PLWNDA(ID,IS,IPTST) = SWINEB*AC2(ID,IS,KCGRD(1))  40.85 40.00
            GENC0(ID,IS,1) = GENC0(ID,IS,1) + SWINEB*AC2(ID,IS,KCGRD(1))  40.85
!
!           *** test output ***
!
            IF (ITEST.GE. 30 .AND. TESTFL) THEN
              WRITE(PRINTF,101) ZARG,ZLOG,ID,IS
  101         FORMAT(' SW4: ZARG ZLOG ID IS  :',2E12.4,' in',2I3)
              WRITE(PRINTF,102) COS1, COS2, XFAC2, BETA
  102         FORMAT(' SW4: COS1/2 XFAC BETA :',4E12.4)
            ENDIF

          ENDIF
        ENDDO
!
       IF (ITEST.GE. 20 .AND. TESTFL) THEN
          WRITE(PRINTF,1102)  SIGMA, WAVEN, CW1, ZCN
 1102     FORMAT(' SW4: SIG WAV CW1 ZCN :',4E12.4)
        ENDIF
      ENDDO
!
!     *** test output ***
!
      IF (ITEST.GE. 20.AND.TESTFL) THEN
        WRITE(PRINTF,9001) KCGRD(1), IDWMIN, IDWMAX
 9001   FORMAT(' SW4: POINT IDWMIN IDWMAX :',3I5)
        WRITE(PRINTF,6053) WIND10,UFRIC,THETAW*180./PI
 6053   FORMAT(' SW4: WIND10 UFRIC THETAW :',3E12.4)
        WRITE(PRINTF,7136) PWIND(9), PWIND(16), PWIND(17)
 7136   FORMAT(' SW4: RHOAW  RHOA  RHOW   :',3E12.4)
        WRITE(PRINTF,7126) PWIND(14), PWIND(15), ZTEN
 7126   FORMAT(' SW4: ALPHA XKAPPA ZTEN   :',3E12.4)

      END IF
!
      RETURN
!     end of subroutine SWIND4
      END
!
!****************************************************************
!
      SUBROUTINE SWIND5 (SPCSIG  ,THETAW  ,ISSTOP  ,
     &                   UFRIC   ,KWAVE   ,IMATRA  ,IDCMIN  ,
     &                   IDCMAX  ,AC2     ,ANYWND  ,PLWNDA  ,
     &                   SPCDIR  ,GENC0                     )
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.41: Marcel Zijlema
!     40.53: Andre van der Westhuysen
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     40.41, Aug. 04: COS(THETA-THETAW) replaced by sumrule to make it cheaper
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.53, Aug. 04: changes parameters of Yan formulae in case of Alves and
!                     Banner whitecapping method
!     40.85, Aug. 08: store wind input for output purposes
!
! 2. Purpose
!
!     Computation of the source term for the wind input for a
!     third generation wind growth model:
!
!     The exponential input term is according to Yan (1987). This
!     input term is valid for the higher frequency part of the
!     spectrum (strongly forced wave components). The expression
!     reduces to the Snyder (1982) expression form for spectral
!     wave components with weak wind forcing and to the Plant (1982)
!     form for more strongly forced wave components:
!
!  3. Method
!
!     The expression reads -->   with  X = Ustar / C
!
!            / /      2                      \
!     Sin = | | 0.04 X + 0.00544 X + 0.000055 | * cos (theta)
!            \ \                             /
!                       \
!              - 0.00031 | sigma * AC2(d,s,x,y)
!                       /
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
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!        FPM         Pierson Moskowitz frequecy (WAM)
!        THETA       Spectral direction
!        THETAW      Mean direction of the relative wind vector
!        UFRIC       Wind friction velocity
!
!        one and more dimensional arrays:
!        ---------------------------------
!        KWAVE     2D    Wavenumber
!        IMATRA    2D    Coefficients of right hand side of matrix
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!        ANYWND    1D    Wind input for bin considered
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SOURCE
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
!        ---
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every spectral bin that fall within a sweep considered
!     compute the source term and store the results in IMATRA
!   --------------------------------------------------------------
!   End of the subroutine SWIND5
!   --------------------------------------------------------------
!
!     10. SOURCE
!
!***********************************************************************
!
      INTEGER  IDDUM  ,ID     ,IS     ,ISSTOP
!
      REAL     UFRIC  ,THETA  ,THETAW ,SIGMA  ,SWINEB ,
     &         CTW    ,STW    ,COSDIF ,                                   40.41
     &         USTAC1 ,USTAC2 ,COF1   ,COF2   ,COF3   ,COF4
!
      REAL    AC2(MDC,MSC,MCGRD)   ,
     &        IMATRA(MDC,MSC)      ,
     &        KWAVE(MSC,ICMAX)     ,
     &        PLWNDA(MDC,MSC,NPTST)
      REAL    GENC0(MDC,MSC,MGENR)                                        40.85
!
      INTEGER IDCMIN(MSC)          ,
     &        IDCMAX(MSC)
!
      LOGICAL  ANYWND(MDC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWIND5')
!
!     *** input according to Yan (1987) ***
!
      COF1 = 0.04
      COF2 = 0.00544
      COF3 = 0.000055
      COF4 = 0.00031
!
!     adapted Yan fit for use of Alves and Banner method                  40.53
!
      IF (IWCAP.EQ.7) THEN                                                40.53
         COF1 = 0.04                                                      40.53
         COF2 = 0.00552                                                   40.53
         COF3 = 0.000052                                                  40.53
         COF4 = 0.000302                                                  40.53
      END IF                                                              40.53
!
      CTW  = COS(THETAW)                                                  40.41
      STW  = SIN(THETAW)                                                  40.41
      DO IS = 1, ISSTOP
        SIGMA  = SPCSIG(IS)
        USTAC1 = ( UFRIC * KWAVE(IS,1) ) / SIGMA
        USTAC2 = USTAC1 * USTAC1
        TEMP3  = ( COF1 * USTAC2 + COF2 * USTAC1 + COF3)
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC, MDC ) + 1
          IF ( ANYWND(ID) ) THEN
            THETA  = SPCDIR(ID,1)                                         30.82
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW                  40.41
            SWINEB = TEMP3 * COSDIF - COF4                                40.41
            SWINEB = MAX ( 0. , SWINEB * SIGMA )
            IMATRA(ID,IS) = IMATRA(ID,IS) + SWINEB * AC2(ID,IS,KCGRD(1))
            IF (TESTFL) PLWNDA(ID,IS,IPTST) = SWINEB*AC2(ID,IS,KCGRD(1))  40.85 40.00
            GENC0(ID,IS,1) = GENC0(ID,IS,1) + SWINEB*AC2(ID,IS,KCGRD(1))  40.85
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST.GE. 60.AND.TESTFL) THEN
        WRITE(PRINTF,6000) KCGRD(1), THETAW*180./PI, UFRIC
 6000   FORMAT(' SWIND5: POINT THETAW UFRIC  :',I5,2E12.4)
        WRITE(PRINTF,6100) COF1, COF2, COF3, COF4
 6100   FORMAT(' SWIND5: COF1 COF2 COF3 COF4 :',4E12.4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     end of subroutine SWIND5
      END
!
