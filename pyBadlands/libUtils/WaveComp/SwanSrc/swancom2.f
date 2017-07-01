!
!     SWAN/COMPU      file 2 of 5
!
!     PROGRAM SWANCOM2.FOR
!
!     This file SWANCOM2 of the main program SWAN
!     include the next subroutines (mainly subroutines for
!     the source terms for dissipation and some general stuff):
!
!     DISSIPATION SOURCE TERMS :
!
!     SBOT    (Bottom friction)
!     SVEG    (Dissipation due to vegetation)                             40.55
!     FRABRE  (Fraction of breaking waves)                                30.77
!     SSURF   (Wave breaking: five formulations)
!     SWCAP   (White capping: seven formulations)                         40.53
!     BRKPAR  (wave breaking criterion according to Nelson (1987))
!     CNTAIL  (contributions to the spectrum of the high frequency tail)
!     PLTSRC  (store the values for plot of the source terms and spec.)
!
!****************************************************************
!
      SUBROUTINE SBOT (ABRBOT  ,DEP2    ,ECOS    ,ESIN    ,AC2     ,      41.04
     &                 IMATDA  ,KWAVE   ,SPCSIG  ,UBOT    ,UX2     ,      30.72
     &                 UY2     ,IDCMIN  ,IDCMAX  ,
     &                 PLBTFR  ,ISSTOP  ,DISSC1  ,VARFR   ,FRCOEF  )      40.67
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
!     20.68: Nico Booij
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!     40.61: Marcel Zijlema
!     40.67: Nico Booij
!     41.04: Marcel Zijlema
!
!  1. Updates
!
!     20.68, Jan. 96: subroutine restructured variable friction coefficient
!                     introduced Putnam model replaced by Collins
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Sep. 06: introduce DISBOT variable for output purposes
!     40.67, Jun. 07: more accurate computation fo dissipation terms
!     41.04, Mar. 09: frequency-dependent JONSWAP formulation
!
!  2. Purpose
!
!     Computation of the source terms due to bottom friction
!
!  3. Method
!
!     In SWAN several bottom friction dissipation models are computed, i.e.:
!
!     IBOT = 1   Jonswap bottom friction model
!     IBOT = 2   Collins bottom friction model
!     IBOT = 3   Madsen bottom friction model (see Tolman)
!
!     Both methods are implemented in SWAN and the user has to make
!     a choice in the input file.
!
!     1. Jonswap model:
!     -----------------
!
!     The bottom interaction term SEbf(s,d) is supposed to take the
!     Jonswap form (Hasselman et al. 1973):
!                         2
!                    sigma  E(s,d)
!     SEbf = -GAMMA ----------------
!                     2     2
!                    g  sinh  (kD)
!                                                                 2 -3
!     where GAMMA is the decay parameter ,(default GAMMA = 0.067 m s  ).
!     In the Jonswap form the current velocities are not taken into
!     account.
!
!     2. COLLINS model:
!     -----------------
!
!     The energy dissipation due to bottom friction is modelled
!     according the quadratic friction law:
!                    2
!      SE = Tau * |U|
!
!      which for a spectrum can be written as:
!                          2
!                     sigma
!      SE(s,d)= - ---------------- * (Cfw.Ub + Cfc.Uc) * E(s,d)
!                       2
!                 g sinh (K(s) * D)
!
!     Ub is the velocity due to the wave at the bottom
!
!     The current velocity is Uc
!
!     2. MADSEN formulation:
!     ----------------------
!
!     The bottom friction dissipation applying Madsen formulation is as
!     follows:
!
!                          fw [n - 1/2] UBR E(s,d)
!     [1]    Sdsb(s,d) = -  ------------------------
!                                      D
!
!     in which :
!                            2
!                           s * D
!     [1a]   (n - 1/2) = -------------
!                                2
!                        2 g sinh (kD)
!
!     UBOT(IX,IY) is computed in the subroutine SINTGRL. The friction
!     factor fw is estimated using the formulation of Jonsson (1963,
!     1966a):
!
!                1                1                        Ab,r
!     [2]     -------- + log  { ---------- } = mf + log  { ----- }
!            4 sqrt(fw)     10  4 sqrt(fw)             10   Kn
!
!     with:
!
!               2        //      1
!     [3]   Ab,r  = 2 * // -------------- E(s,d) ds dd
!                      //      2
!                          sinh (kD)
!
!     with: Ab,r is the representative near bottom excursion
!                amplitude
!           Kn   equivalent roughness
!           mf   constant ( mf = -0.08) (determined by Jonsson
!                                        and Carlssen 1976 )
!
!     [2] is only valid for Ab,r/Kn larger than approximately 1.
!     For smaller values a constant value of fw is used (fw = 0.3
!     for Ab,r/Kn < 1.57 )
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     INTEGERS :
!     --------
!
!     IX          Counter of gridpoints in x-direction
!     IY          Counter of gridpoints in y-direction
!     IS          Counter of relative frequency band
!     ID          Counter of the spectral direction
!     IBOT        Indicator if bottom friction is on
!     ICUR        Indicator if a current is present
!     ITER        Number of iteration i.e. number of full sweeps
!     MBOT        Maximum array size for the array PBOT
!     MXC         Maximum counter of gridppoints in x-direction
!     MYC         Maximum counter of gridppoints in y-direction
!     MSC         Maximum counter of relative frequency
!     MDC         Maximum counter of directional distribution
!     ISSTOP      Maximum counter of wave component in frequency
!                 space that is propagated
!
!     REALS:
!     ---------
!
!     ABRBOT      Near bottom excursion amplitude
!     FACB        an auxiliary factor contributing to bottom friction
!     FW          Friction factor
!     GRAV        Gravitational acceleration
!     KD          Wavenumber * Depth
!     SBOTEO      Sourceterm for the bottom friction to be stored
!                 in the array IMATDA
!     CURR        Main current velocity
!     UC          Absolute value of the current
!
!     one and more dimensional arrays:
!     ---------------------------------
!
!     AC2       2D    Action density
!     DEP2      2D    Depth
!     ESIN      1D    Sin per spectral direction (id)
!     ECOS      1D    Cos per spectral direction (id)
!     IMATDA    2D    Coefficients of diagonal of matrix
!     KWAVE     2D    Wavenumber function of frequency and IC
!     PBOT      1D    Coefficient for bottom friction models
!     UBOT      2D    Near bottom velocity as function of X,Y
!     UX2       2D    Current velocity in y direction as function of X,Y
!     UY2       2D    Current velocity in y direction as function of X,Y
!     DISSC1    2D    Dissipation coefficient, function of sigma and theta
!     FRCOEF    2D    Spatially variable friction coefficient             20.68
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE
!
! 10. Error Messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ------------------------------------------------------------
!     Compute CFBOT according to friction model
!     For every spectral frequency do
!         compute SBOTEO = CFBOT * (sigma/sinh(kd))**2
!         For every spectral direction do
!             add SBOTEO to matrix (IMATDA)
!     -------------------------------------------------------------
!
! 13. Source text
!
      INTEGER  ID     ,IS     ,ISSTOP
!
      REAL     XDUM   ,KD     ,SBOTEO ,FACB  ,
     &         CFW    ,FW     ,CURR   ,UC    ,ABRBOT,
     &         ADUM   ,CDUM   ,DDUM
      REAL     CFBOT(MSC)
      REAL     DSP    ,ETOT   ,EEX    ,EEY   ,EAD
!
      LOGICAL  VARFR
!
      REAL     AC2(MDC,MSC,MCGRD)        ,                                41.04
     &         DEP2(MCGRD)               ,
     &         ECOS(MDC)                 ,
     &         ESIN(MDC)                 ,
     &         IMATDA(MDC,MSC)           ,
     &         KWAVE(MSC,ICMAX)          ,
     &         PLBTFR(MDC,MSC,NPTST)     ,                                40.00
     &         UBOT(MCGRD)               ,
     &         UX2(MCGRD)                ,
     &         UY2(MCGRD)                ,
     &         DISSC1(MDC,MSC,1:MDISP)   ,                                40.67
     &         FRCOEF(MCGRD)                                              20.68
!
      INTEGER  IDCMIN(MSC)               ,
     &         IDCMAX(MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SBOT')
!
      IF ( IBOT .GE. 1 .AND. DEP2(KCGRD(1)) .GT. 0.) THEN
        IF (IBOT.EQ.1) THEN
!
!         *** Jonswap model ***
!
!         PBOT(3) = GAMMA (a) in the Jonswap equation,
!
          CFBOT = PBOT(3) / GRAV**2
        ELSE IF (IBOT.EQ.2) THEN
!
!         *** Collins model ***
!
!         PBOT(2) = [cfw]
!
          IF (VARFR) THEN                                                 20.68
            CFW = FRCOEF(KCGRD(1))
          ELSE
            CFW = PBOT(2)
          ENDIF
          CFBOT = CFW * UBOT(KCGRD(1)) / GRAV
        ELSE IF (IBOT.EQ.3) THEN
!
!             *** Madsen model ***
!
          IF (VARFR) THEN                                                 20.68
            AKN = FRCOEF(KCGRD(1))
          ELSE
            AKN = PBOT(5)
          ENDIF
!
!           *** PBOT(4) = Mf                      ***
!           *** AKN = PBOT(5) = [kn]  (roughness) ***
!
          IF ( (ABRBOT / AKN ) .GT. 1.57 ) THEN
            XDUM = PBOT(4) + LOG10 ( ABRBOT / AKN )
!
!               *** solving the implicit equation using a Newton ***
!               *** Rapshon iteration proces : a + log a = b     ***
!               *** the start value for ADUM = 0.3 because 0.3626 ***
!               *** is the minimum value of ADUM with b=-0.08.    ***
!
            ADUM = 0.3
            DO 28 J = 1, 50
              CDUM  = ADUM
              DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) /
     &                                          ( 1.+ ( 1. / ADUM) )
              ADUM  = ADUM - DDUM
              IF ( ABS(CDUM - ADUM) .LT. 1.E-4 ) GOTO 29
  28        CONTINUE
            WRITE(*,*) ' error in iteration fw: Madsen formulation'
  29        CONTINUE
!                                                 1               1
!               *** computation of FW -->  A = ----- --> FW = -----
!                                              4 uFW          16 A**2
            FW = 1. / (16. * ADUM**2)
          ELSE
            FW = 0.3
          ENDIF
          CFBOT =  UBOT(KCGRD(1)) * FW / (SQRT(2.) * GRAV)
        ELSE IF ( IBOT.EQ.4 ) THEN
!
!            *** Jonswap model with variable friction coefficient  ***
!            *** as function of frequency-dependent directional    ***
!                spreading (varies linearly between 0.038 - 0.067) ***
!
          DO IS = 1, MSC
             ETOT = 0.
             EEX  = 0.
             EEY  = 0.
             DO ID = 1, MDC
                EAD  = SPCSIG(IS)*AC2(ID,IS,KCGRD(1))
                ETOT = ETOT + EAD
                EEX  = EEX  + EAD * ECOS(ID)
                EEY  = EEY  + EAD * ESIN(ID)
             ENDDO
             IF ( ETOT.GT.0. ) THEN
                XDUM = 1.-MIN(1.,SQRT(EEX*EEX+EEY*EEY)/ETOT)
                DSP  = SQRT(2.*XDUM) *180./PI
             ELSE
                DSP  = 0.
             ENDIF
             IF ( DSP.LT.PBOT(8) ) THEN
                CFBOT(IS) = PBOT(6)
             ELSEIF ( DSP.GT.PBOT(9) ) THEN
                CFBOT(IS) = PBOT(7)
             ELSE
                CFBOT(IS) = PBOT(6) + (PBOT(7)-PBOT(6))*(DSP-PBOT(8))/
     &                                (PBOT(9)-PBOT(8))
             ENDIF
             CFBOT(IS) = CFBOT(IS) / GRAV**2
          ENDDO
        ENDIF
!
!       *** test output ***
!
        IF (TESTFL .AND. ITEST.GE.60) THEN
          WRITE (PRTEST, 910) IBOT, KCGRD(1), DEP2(KCGRD(1)), CFBOT(1)
 910      FORMAT (' SBOT :IBOT INDX DEP CFBOT:', 2I5, 2E12.4)
        END IF
!
        DO 700 IS = 1, ISSTOP
          KD = KWAVE(IS,1) * DEP2(KCGRD(1))
          IF ( KD .LT. 10. ) THEN
            FACB = CFBOT(IS) * (SPCSIG(IS) / SINH(KD)) **2                41.04 40.57 30.72
!
            DO 690 IDDUM = IDCMIN(IS) , IDCMAX(IS)
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
              SBOTEO = FACB                                               40.57
              IF (IBOT.EQ.2 .AND. ICUR.EQ.1 .AND. PBOT(1).GT.0.) THEN
!               additional dissipation due to current, seldom used
                CURR = UX2(KCGRD(1))*ECOS(ID) + UY2(KCGRD(1))*ESIN(ID)
                UC   = ABS(CURR)
!               PBOT(1) = [cfc]
                SBOTEO = FACB + PBOT(1) * UC *                            40.57
     &                         (SPCSIG(IS) / SINH(KD)) **2                30.72
              END IF
!
!             *** store the results in the array IMATDA             ***
!             *** if testfl store results in array for isoline plot ***
!
              IMATDA(ID,IS) = IMATDA(ID,IS) + SBOTEO
              IF (TESTFL) PLBTFR(ID,IS,IPTST) = -1.* SBOTEO               40.00
              DISSC1(ID,IS,3) = DISSC1(ID,IS,3) + SBOTEO                  40.67
 690        CONTINUE
          END IF
 700    CONTINUE
!
      ENDIF
!
!     End of subroutine SBOT
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE SVEG ( DEP2   ,IMATDA   ,ETOT   ,SMEBRK    ,
     &                  KMESPC ,PLVEGT   ,
     &                  IDCMIN ,IDCMAX   ,ISSTOP ,DISSC1    ,
     &                  NPLA2  )
!
!****************************************************************
!
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4
      USE M_GENARR
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
!     40.55: Bastiaan Burger, Martijn Meijer
!     40.61: Marcel Zijlema
!     40.58: Tomo Suzuki, Marcel Zijlema
!
!  1. Updates
!
!     40.55, May. 05: implementation of vegetation dissipation formula
!     40.61, Sep. 06: introduce DISVEG variable for output purposes
!     40.58, Nov. 08: some modifications and corrections
!
!  2. Purpose
!
!     Computation of the source term due to vegetation dissipation
!
!  3. Method
!
!     The energy dissipation due to vegetation is described by a
!     Morrison type equation, modelling the plants as vertical,
!     noncompliant cylinders, neclecting swaying motions induced by
!     waves. Vegetation characteristics that are used as input are
!     drag coefficient, vegetation height, plant density and diameter.
!
!     The formula used in SWAN is due to Dalrymple (1984)
!     (see Mendez and Losada, 2004):
!
!     d(Cg E)
!     ------- = -epsv
!        dx
!
!     with dissipation due to vegetation:
!
!     epsv = 1/2 * 1/sqrt(pi) * rho * Cd * bv * Nv *
!                                                                3
!                              (gk/2sigma)^3 * ((A + B)/C) * Hrms
!
!     where
!
!     rho   = water density
!     k     = wave number
!     sigma = angular frequency
!     Cd    = drag coefficient
!     bv    = stem thickness
!     Nv    = vegetation density
!     Hrms  = rms wave height
!
!     and the coefficients
!
!     A = sinh^3 k*ah
!     B = 3*sinh k*ah
!     C = 3k*cosh^3 kh
!
!     with ah the vegetation height
!
!     The source term to be used in SWAN is based on the
!     corresponding dissipation rate and reads
!
!     Dtot = -epsv / rho / g
!
!     In the formulation, the mean average wavenumber according
!     to the WAM-formulation and the mean frequency will be employed
!
!     Now, the source term is:
!
!                      E
!     Sveg = Dtot *  ------ = factor * sqrt(Etot) * E
!                     Etot
!
!     and is linearized by means of the Picard iteration
!
!
!  4. Argument variables
!
!     DEP2        water depth
!     DISSC1      dissipation coefficient
!     ETOT        total energy per spatial gridpoint
!     IDCMIN      frequency dependent counter in directional space
!     IDCMAX      frequency dependent counter in directional space
!     IMATDA      coefficients of diagonal of matrix
!     ISSTOP      maximum counter of wave component in frequency
!                 space that is propagated
!     KMESPC      mean average wavenumber according to the WAM-formulation
!     NPLA2       number of plants per square meter (depth-averaged)
!     PLVEGT      array containing the vegetation source term for test-output
!     SMEBRK      mean frequency according to first order moment
!
      INTEGER ISSTOP, IDCMIN(MSC), IDCMAX(MSC)
      REAL    DEP2(MCGRD)          ,
     &        IMATDA(MDC,MSC)      ,
     &        DISSC1(MDC,MSC,MDISP),
     &        PLVEGT(MDC,MSC,NPTST),
     &        NPLA2 (MCGRD)
      REAL    ETOT, SMEBRK, KMESPC
!
!  6. Local variables
!
!     A     :     auxiliary variable
!     B     :     auxiliary variable
!     C     :     auxiliary variable
!     D     :     auxiliary variable
!     ID    :     counter of the spectral direction
!     IDDUM :     counter
!     IENT  :     number of entries
!     IK    :     counter
!     IL    :     counter
!     IS    :     counter of relative frequency band
!     KD    :     wavenumber times water depth
!     KVEGH :     wavenumber times plant height
!     LAYPRT:     part of layer below water level
!     SINHK :     sinh(kh)
!     SLAYH :     total sum of layer thicknesses
!     SLAYH1:     sum of layer thicknesses below water level
!     SLAYH2:     sum of layer thicknesses below water level
!     SVEG1 :     layer-independent dissipation factor
!     SVEG2 :     total sum of dissipation factor over layers
!     SVEGET:     source term containing dissipation due to vegetation
!                 to be stored in the array IMATDA
!
      INTEGER ID, IDDUM, IENT, IK, IL, IS
      REAL    A, B, C, D, KD, KVEGH, LAYPRT, SINHK, SLAYH,
     &        SLAYH1, SLAYH2, SVEG1, SVEG2, SVEGET
!
!  9. Subroutines calling
!
!     SOURCE
!
! 12. Structure
!
!     Vegetation parameters are given per layer thickness, so for
!     each layer the contribution to wave damping is calculated
!
!     This routine checks in which layer the water level is present
!
!            -----------
!                        ILMAX
!            -----------
!              _ d       2
!            --|--------
!              |         1
!            --|--------
!
!     d     = waterdepth
!     ILMAX = number of layers in grid point
!
!     Subsequently, the vegetation parameters up to the layer where the
!     water level is in, are used to calculate dissipation for each layer
!
!     Thereafter, the contributions to disspation are summed up
!
!     With this summation the total dissipation due to vertical varying
!     vegetation is calculated
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SVEG')

!     --- compute layer-independent vegetation dissipation factor

      KD    = KMESPC * DEP2(KCGRD(1))
      IF ( KD.GT.10. ) RETURN
      C     = 3.*KMESPC*(COSH(KD))**3
      SVEG1 = SQRT(2./PI) * GRAV**2 * (KMESPC/SMEBRK)**3 * SQRT(ETOT)
     &                    * NPLA2(KCGRD(1)) / C

!     --- compute dissipation factor for each layer and summed up

      SLAYH = 0.
      DO IL = 1, ILMAX
         SLAYH = SLAYH + LAYH(IL)
      ENDDO

      KVEGH = 0.
      C     = 0.
      D     = 0.
      SVEG2 = 0.

      IF ( DEP2(KCGRD(1)).GT.SLAYH ) THEN

         DO IL = 1, ILMAX
            KVEGH = KVEGH + KMESPC * LAYH(IL)
            SINHK = SINH(KVEGH)
            A     = C
            B     = D
            C     = SINHK**3
            D     = 3.*SINHK
            A     = C - A
            B     = D - B
            SVEG2 = SVEG2 + VEGDRL(IL)*VEGDIL(IL)*VEGNSL(IL)*(A + B)
         END DO

      ELSE IF ( DEP2(KCGRD(1)).LT.LAYH(1) ) THEN

         SINHK = SINH(KD)
         A     = SINHK**3
         B     = 3.*SINHK
         SVEG2 = VEGDRL(1)*VEGDIL(1)*VEGNSL(1)*(A + B)

      ELSE

         SLAYH1 = 0.
         SLAYH2 = 0.
         LAYPRT = 0.
         VGLOOP : DO IL = 1, ILMAX
            SLAYH1 = SLAYH1 + LAYH(IL)
            IF (DEP2(KCGRD(1)).LE.SLAYH1) THEN
               DO IK = 1, IL-1
                  SLAYH2 = SLAYH2 + LAYH(IK)
               END DO
               LAYPRT = DEP2(KCGRD(1)) - SLAYH2
               DO IK = 1, IL-1
                  KVEGH = KVEGH + KMESPC * LAYH(IK)
                  SINHK = SINH(KVEGH)
                  A     = C
                  B     = D
                  C     = SINHK**3
                  D     = 3.*SINHK
                  A     = C - A
                  B     = D - B
                  SVEG2 = SVEG2+VEGDRL(IK)*VEGDIL(IK)*VEGNSL(IK)*(A + B)
               END DO
               KVEGH = KVEGH + KMESPC * LAYPRT
               SINHK = SINH(KVEGH)
               A     = C
               B     = D
               C     = SINHK**3
               D     = 3.*SINHK
               A     = C - A
               B     = D - B
               SVEG2 = SVEG2 + VEGDRL(IL)*VEGDIL(IL)*VEGNSL(IL)*(A + B)
               EXIT VGLOOP
            END IF
         END DO VGLOOP

      END IF

!     --- compute total dissipation

      SVEGET = SVEG1 * SVEG2
!
!     *** test output ***
!
      IF (TESTFL .AND. ITEST.GE.60) THEN
         WRITE (PRTEST, 110) IVEG, KCGRD(1), DEP2(KCGRD(1)), SVEGET
 110     FORMAT (' SVEG :IVEG INDX DEP VEGFAC:', 2I5, 2E12.4)
      END IF

      DO IS = 1, ISSTOP
         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1

!           *** store the results in the array IMATDA ***
!           *** if testfl store results in array for isoline plot ***

            IMATDA(ID,IS) = IMATDA(ID,IS) + SVEGET
            IF (TESTFL) PLVEGT(ID,IS,IPTST) = -1.* SVEGET
            DISSC1(ID,IS,5) = DISSC1(ID,IS,5) + SVEGET

         END DO
      END DO

      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE FRABRE ( HM, ETOT, QBLOC )                               30.77
!
!****************************************************************
!
      USE SWCOMM4                                                         40.41
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
!     30.77: Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.77, Sep. 98: the discontinuity at B = 0.9 has been removed and
!                     the discontinuity at B = 0.3 is changed in a discontinuity
!                     at B = 0.2 for which QBLOC = 1.E-9
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     to compute the fraction of breaking waves in point ix,iy
!     of the computational grid
!
!  3. Method (updated...)
!
!      The fraction of breaking waves in a point ix,iy is given by
!      the implicit relation:
!
!        1 - Qb        ETOT
!        ------ = -8 * -----
!        ln Qb         HM**2
!
!        from which Qb can be found by solving the equation:
!
!                         ETOT
!        F = 1 - Qb + 8 * ----  * ln(Qb) = 0.
!                           2
!                         HM
!
!        The following appproximation is applied:
!
!                            2
!  (1)|   B = sqrt( 8 ETOT/HM ), i.e. B = Hrms/HM
!
!
!     |   Qo = 0.                                      B <= 0.5
!  (2)|                 2
!     |   Qo = ( 2B -1 )                         0.5 < B <= 1
!
!
!     applying the Newton-Raphson procedure (for 0.2<B<1.0):
!
!     |   Qb = 0.                                      B <= 0.2
!     |
!     |                                   2
!     |               2  Qo - exp((Qo-1)/B )
!  (3)|   Qb = Qo  - B   ------------------      0.2 < B <  1.0
!     |                   2               2
!     |                  B  - exp((Qo-1)/B )
!     |
!     |
!     |   Qb = 1.                                      B >= 1.0
!     |
!
!     Here the parameters ETOT and HM are determined in the subroutine
!     SINTGRL
!
!  4. Argument variables
!
!     ETOT    input  total energy per spatioal gridpoint
!     HM      input  maximum wave height
!     QBLOC   output second iteration of the fraction of breaking waves
!
      REAL    ETOT,  HM,  QBLOC
!
!  5. Parameter variables
!
!  6. Local variables
!
!     B       dummy variable
!     B2      dummy variable: B**2
!     IENT    number of entries
!     QO      first estimate of the fraction of breaking waves
!     Z       dummy variable
!
      INTEGER IENT
      REAL    B,  B2,  QO,  Z
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SINTGRL
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!   ------------------------------------------------------------
!   Read the total wave energy ETOT and the maximum waveheight HM
!     If HM > 0. and ETOT > 0., then
!       Compute factor B according to equation (1)
!     Else
!       B = 0
!     ------------------------------------------------------------
!     Compute first estimate Qo according to equation (2)
!     Compute Qb according to equation (3)
!   ------------------------------------------------------------
!   End of FRABRE
!   ------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'FRABRE')
!
      IF ( (HM .GT. 0.) .AND. (ETOT .GE. 0.) ) THEN
        B = SQRT(8. * ETOT / (HM*HM) )
      ELSE
        B = 0.0
      END IF
!
      IF ( B .LE. 0.5 ) THEN
        QO = 0.
      ELSE IF ( B .LE. 1.0 ) THEN
        QO = (2.*B - 1.)**2
      END IF
!
      IF ( B .LE. 0.2 ) THEN
        QBLOC = 0.0
      ELSE IF ( B .LT. 1.0 ) THEN
!
!       *** second iteration to find Qb ***
!
        B2 = B*B
        Z  = EXP((QO-1.)/B2)
        QBLOC = QO - B2 * (QO-Z)/(B2-Z)
      ELSE
        QBLOC = 1.0
      END IF
!
      IF ( TESTFL .AND. ITEST .GE. 110 ) THEN
        WRITE (PRINTF,6120) ETOT, HM, B, QBLOC
 6120   FORMAT (' FRABRE: ETOT  HM  B  QB     : ',4E12.4)
      END IF
!
!     End of subroutine FRABRE
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE SSURF (ETOT    ,HM      ,QB      ,SMEBRK  ,              30.81
     &                  KMESPC  ,SPCSIG  ,AC2     ,IMATRA  ,              30.81
     &                  IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,              30.81
     &                  ISSTOP  ,DISSC0  ,DISSC1  )                       40.67 40.61 30.81 30.21
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
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
!     30.62: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!     40.61: Marcel Zijlema
!     40.67: Nico Booij
!     41.03: Andre van der Westhuysen
!     41.06: Gerbrant van Vledder
!
!  1. Updates
!
!     30.62, Aug. 97: Prevented a possible division by zero
!     30.82, Sep. 98: Changed indices of PLWBRK-array declaration
!     30.82, Oct. 98: Made subroutine intrinsic DOUBLE PRECISION
!     30.81, Sep. 99: Argumentlist reduced
!     40.13, Jan. 01: PLWBRK corrected (dissipation test output)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Sep. 06: introduce DISSRF variable for output purposes
!     40.67, Jun. 07: more accurate computation of dissipation terms
!     41.03, Feb. 09: extension to alternative surf breaking formula's
!     41.06, Mar. 09: extension to frequency dependent surf breaking
!
!  2. Purpose
!
!     Computation of the source term due to wave breaking with one of the
!     following formulation:
!                 1) Battjes and Janssen (1978)
!                 2) Thornton and Guza (1983)
!
!     Note: white capping is not taken into account
!
!  3. Method
!
!     Basically, the source term for surf breaking is implemented following
!     the approach of Battjes/Janssen (1978) for the energy dissipation:
!
!             Alpha      -     2                  -   SMEBRK
!     Dtot =  ----  Qb * f * Hm              with f = ------
!              4                                      2 * Pi
!
!     Now, the source term is:
!
!                      SIGMA * AC2(ID,IS,IX,IY)
!     Sbr =    Dtot *  ------------------------  =
!                              Etot
!
!
!              Alpha * SMEBRK * Qb * Hm * Hm    SIGMA * AC2(ID,IS,IX,IY)
!         =    ------------------------------ * -------------------------
!                       8 * Pi                            Etot
!
!
!         =  WS * SIGMA * AC2(ID,IS,IX,IY)   =  WS * E
!
!
!     with
!
!     Alpha = PSURF(1)                            ;
!
!                   SMEBRK Qb
!     WS    = Alpha ------ --                     ;
!                     Pi   BB
!                        2
!     BB    = 8 Etot / Hm  = - (1 - Qb) / ln (Qb) ;
!
!
!     The local maximum wave height Hm and mean frequency SMEBRK are computed
!     in subroutine SINTGRL.
!     The fraction of breaking waves Qb is calculated in the subroutine FRABRE
!
!     The new value for the dissipation is computed implicitly using
!     the last computed value for the action density Nold (at the spatial
!     gridpoint under consideration).
!
!     Sbr = WS * N
!
!         = Sbr_new + (d Sbr/d N) (Nnew - Nold)
!
!         = WS * Nnew + SbrD * (Nnew - Nold)
!
!         = (WS + SbrD)* Nnew - SbrD * Nold
!
!         = SURFA1 * Nnew - SURFA0 * Nold
!
!     In order to do this we need the derivative
!     of the source term Sbr to the action density N
!
!             d Sbr     d WS
!     SbrD =  -----  =  ---- * N + WS
!             d N       d N
!
!     Since BB and N are proportional, we have
!
!     d Sbr     d WS                   SMEBRK  (d Qb/ d BB) *BB - Qb
!     -----  =  ---- * BB + WS = Alpha ------  --------------------- * BB + WS
!     d N       d BB                     Pi           sqr(BB)
!
!
!                    SMEBRK d Qb
!            = Alpha ------ ----
!                     Pi    d BB
!
!     With:
!
!     d Qb         1
!     ---- = -------------                 ;
!     d BB   (d BB / d Qb)
!
!                      2
!     d Qb           ln (Qb)
!     ---- = ---------------------------
!     d BB   ln (Qb) + (1 - Qb) (1 / Qb)
!
!            Qb (1 - Qb)
!          = ------------                  ;
!            BB (BB - Qb)
!
!     Hence,
!
!     d Sbr       1 - Qb
!     ----- = WS  -------
!     d N         BB - Qb
!
!
!     Alternatively, the source term for surf breaking is implemented following
!     the approach of Thornton and Guza (1983) for energy dissipation:
!
!               3  -
!              B * f                3                          -   SMEBRK
!     Dtot =  ------- * INT(0,inf){H * W(H) * p(H)}dH     with f = ------
!              4 * d                                               2 * Pi
!
!                    3
!               3 * B * SMEBRK             3
!          =  ------------------ * W * Hrms
!              32 * sqrt(Pi) * d
!
!                           3                             3
!          with INT(0,inf){H * p(H)}dH = 3/4*sqrt(Pi)*Hrms
!
!
!          after Thornton and Guza (1983), Eqs. (24),(25)
!
!     and
!
!           Hrms  n
!     W = (------)
!           Hmax
!
!
!     For implementation details, see the Scientific/Technical documentation.
!
!
!  4. Argument variables
!
!     AC2     input :   Action density array
!     DISSC0  output:   Dissipation coefficient as explicit part
!                       (meant for output)
!     DISSC1  output:   Dissipation coefficient as implicit part
!                       (meant for output)
!     ETOT    input :   Total energy per spatial gridpoint
!     HM      input :   Maximum wave height
!     ICMAX   input :   Maximum number of elements in KCGRD array
!     IDCMIN  input :   Minimum number for counter IDDUM
!     IDCMAX  input :   Maximum number for counter IDDUM
!     IMATDA  output:   Coefficient of diagonal matrix (2D)
!     IMATRA  output:   Coefficient of righthandside of matrix
!     ISSTOP  input :   Maximum for counter IS
!     KMESPC  input :   Mean average wavenumber according to the WAM-formulation
!     PLWBRK  output:   array containing the surf breaking source term
!                       for test-output
!     QB      input :   Fraction of breaking waves
!     SMEBRK  input :   Mean frequency according to first order moment
!
      INTEGER        ISSTOP,
     &         IDCMIN(MSC), IDCMAX(MSC)
!
      REAL     AC2(MDC,MSC,MCGRD)   ,
     &         DISSC0(MDC,MSC,MDISP),                                     40.67
     &         DISSC1(MDC,MSC,MDISP),                                     40.67
     &         IMATDA(MDC,MSC)      ,
     &         IMATRA(MDC,MSC)      ,
     &         PLWBRK(MDC,MSC,NPTST)                                      40.00
      REAL     SPCSIG(MSC)
!
      REAL     ETOT,  HM,  QB, SMEBRK, KMESPC                             30.81
!
!  5. Parameter variables
!
!  6. Local variables
!     BB      Rate between the total energy and the energy
!             according to the maximum wave height HM
!     DIS0    Dummy variable
!     ID      Counter for directional steps
!     IDDUM   Counter
!     IENT    Number of entries
!     IS      Counter for frequency steps
!     SURFA0  Coefficient for old source term in matrix equation
!             (i.e. SURFA0 * Nold = right hand side of matrix equation)
!     SURFA1  Coefficient for new source term in matrix equation
!     WS      Wavebreaking source term coefficient = DTOT/ETOT
!     SbrD    Derivative of source term for surf breaking (Sbr) to action density
!
      INTEGER          ID,       IDDUM,   IENT,   IS
      REAL             PP,       FAC,     EPTOT,  ETOT0,
     &                 ECS(MDC), FMIN,    FMAX,   FRFAC(MSC)
      DOUBLE PRECISION BB,       DIS0,    SbrD,
     &                 SURFA0,   SURFA1,  WS  ,
     &                 TEMP1 ,   TEMP2
      REAL             SwanIntgratSpc
!
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE
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
!     ------------------------------------------------------------
!     Get HM, QB and ETOT from the subroutine SINTGRL
!     For spectral direction IS and ID do,
!       get the mean energy frequency average over the full spectrum
!       If ETOT > 0 then
!         compute source term for energy dissipation SURFA0 and SURFA1
!       Else
!         source term for wave breaking is 0.
!       End if
!       ----------------------------------------------------------
!       Compute source terms for energy averaged frequency
!       Store results in the arrays IMATDA and IMATRA
!     ------------------------------------------------------------
!     End of SSURF
!     -------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SSURF')
!
!     ALFA = PSURF(1)   <default = 1.0>
!
      BB = 8D0 * DBLE(ETOT) / ( DBLE(HM)**2 )                             41.03 30.82
      SURFA0 = 0D0
      SURFA1 = 0D0
!
      IF ( ISURF.LE.3 ) THEN                                              41.03
!
!        --- Battjes and Janssen (1978)
!
         IF (REAL(BB) .GT. 0. .AND.                                       30.82
     &       REAL(ABS(BB - DBLE(QB))) .GT. 0.) THEN                       30.82
            IF ( BB .LT. 1D0 ) THEN                                       41.03
               WS  = ( DBLE(PSURF(1)) / DBLE(PI)) *                       30.82
     &                 DBLE(QB) * DBLE(SMEBRK) / BB                       30.82
               SbrD = WS * (1D0 - DBLE(QB)) / (BB - DBLE(QB))             41.03 30.82 40.00
            ELSE
               WS  = ( DBLE(PSURF(1)) / DBLE(PI)) * DBLE(SMEBRK)          30.82
               SbrD = 0D0
            END IF
            SURFA0 = SbrD
            SURFA1 = WS + SbrD
         ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
         ENDIF
!
      ELSEIF ( ISURF.EQ.4 ) THEN                                          41.03
!
!        --- Thornton and Guza (1983)
!
         IF ( BB.GT.0D0 ) THEN
            IF ( BB.LT.1D0 ) THEN
               WS = 75D-2*DBLE(PSURF(2))*DBLE(PSURF(1))**3*DBLE(SMEBRK)*
     &              BB**(0.5*(PSURF(4)+1))/DBLE(SQRT(PI))
               SbrD = 5D-1*DBLE(3.+PSURF(4))*WS
            ELSE
               WS = 75D-2*DBLE(PSURF(2))*DBLE(PSURF(1))**3*DBLE(SMEBRK)/
     &              DBLE(SQRT(PI))
               SbrD = WS
            ENDIF
            SURFA0 = SbrD - WS
            SURFA1 = SbrD
         ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
         ENDIF
!
      ENDIF
!
!     *** store the results for surf wave breaking  ***
!     *** in the matrices IMATDA and IMATRA         ***
!
      FRFAC = 1.                                                          41.06
      IF (IFRSRF.EQ.1) THEN                                               41.06
         PP    = PSURF(11)
         FMIN  = PI2*PSURF(12)
         FMAX  = PI2*PSURF(13)
         ECS   = 1.
         ETOT0 = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, ECS, SPCSIG,
     &                          ECS, 0., 0., AC2(1,1,KCGRD(1)), 1 )
         EPTOT = SwanIntgratSpc(PP, FMIN, FMAX, SPCSIG, ECS, SPCSIG,
     &                          ECS, 0., 0., AC2(1,1,KCGRD(1)), 1 )
         FAC   = ETOT0/EPTOT
         IF ( ETOT0.GT.1.E-8 ) THEN
            DO IS = 1, ISSTOP
               FRFAC(IS) = FAC*SPCSIG(IS)**PP
            END DO
         END IF
      END IF
      TEMP1 = SURFA0
      TEMP2 = SURFA1
      DO 101 IS = 1, ISSTOP
        SURFA0 = TEMP1*FRFAC(IS)                                          41.06
        SURFA1 = TEMP2*FRFAC(IS)                                          41.06
        DO 100 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IMATDA(ID,IS) = IMATDA(ID,IS) + REAL(SURFA1)                    30.82
          DIS0 = SURFA0 * DBLE(AC2(ID,IS,KCGRD(1)))                       30.82
          IMATRA(ID,IS) = IMATRA(ID,IS) + REAL(DIS0)                      30.82
          IF (TESTFL) PLWBRK(ID,IS,IPTST) = REAL(SURFA0-SURFA1)           40.13
          DISSC0(ID,IS,2) = DISSC0(ID,IS,2) - REAL(DIS0)                  40.67 30.82
          DISSC1(ID,IS,2) = DISSC1(ID,IS,2) + REAL(SURFA1)                40.67 30.82
 100    CONTINUE
 101  CONTINUE
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 110 ) THEN
        WRITE(PRINTF,6021) SURFA1,SURFA0
 6021   FORMAT (' SSURF : SURFA1 SURFA0     :',2D12.4)
        WRITE(PRINTF,6020) HM, QB, ETOT, SMEBRK
 6020   FORMAT ('       : HM QB ETOT SMEBRK :',4E12.4)
      END IF
!
!
!     end of the subroutine SSURF
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE SWCAP  (SPCDIR  ,SPCSIG  ,KWAVE   ,AC2     ,             40.02
     &                   IDCMIN  ,IDCMAX  ,ISSTOP  ,                      40.02
     &                   ETOT    ,IMATDA  ,IMATRA  ,PLWCAP  ,             40.02
     &                   CGO     ,UFRIC   ,                               40.53
     &                   DEP2    ,DISSC1  ,DISSC0  )                      40.67 40.61 40.12
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_WCAP
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
!     40.02: IJsbrand Haagsma
!     40.12: IJsbrand Haagsma
!     40.30: Gerbrant van Vledder
!     40.31: Gerbrant van Vledder
!     40.41: Gerbrant van Vledder
!     40.41: Marcel Zijlema
!     40.53: Andre van der Westhuysen
!     40.61: Marcel Zijlema
!     40.63: Andre van der Westhuysen
!     40.67: Nico Booij
!
!  1. Updates
!
!     40.02, Jan. 00: New, based on the old SWCAP1-5 subroutines
!     40.12, Nov. 00: Added WCAP to dissipation output (bug fix)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.53: Aug. 04: white-capping based Alves and Banner (2003) method
!     40.61, Sep. 06: introduce DISWCP variable for output purposes
!     40.63, Apr. 07: a correction to Alves and Banner method
!     40.67, Jun. 07: more accurate computation of source terms
!                     DISSIP and DISIMP renamed DISSC1 and DISSC0 (as elsewhere)
!
!
!  2. Purpose
!
!     Calculates the dissipation due to whitecapping
!
!  3. Method
!
!     Whitecapping dissipation is formulated as follows:
!
!     S_wc(sig,th) = - C_wc E(sig,th)
!
!     where the coefficient C_wc has four basic forms:
!
!     C_wc1 = C_K  sig~ (k/k~): According to Komen (generalised)
!     C_wc2 = C_BJ sig~ (k/k~): According to Battjes-Janssen (modified)
!     C_wc3 = C_LH            : According to Longuet Higgins (1969)
!
!     In these formulations C_K is defined as (Komen; 1994 p. 145):
!
!                                       n1           n2
!     C_K = C1 [(1-delta) + delta (k/k~)  ] (S~/S_PM)
!
!     where C1, delta, n1 and n2 can be varied
!
!
!     C_BJ is defined as:
!
!                    2
!            alpha Hm Qb
!     C_BJ = -----------
!              8 Pi m0
!
!     where alpha can be varied.
!
!     for Hrms > Hm the formulation changes in a limit to (Hrms->Hm; Qb->1):
!
!            alpha
!     C_BJ = -----
!             Pi
!
!     and C_LH is defined as (Komen; 1994 p. 151-152):
!
!                            4   2                        2
!     C_LH = C3 Sqrt[(m0 sig0 )/g ] exp(A) sig0 (sig/sig0)
!
!     where
!                    2    2         4                2       2
!     A = -1/8 [1-eps ] [g /(m0 sig0 )]  with  [1-eps ] = [m2 ] / [m0 m4]
!
!     and C3 can be varied
!
!
!     In these equations the variables have the following meaning:
!
!     Hm   : Maximum wave height
!     Hrms : Root mean square of the wave heights
!     eps^2: Measure for the spectral bandwidth
!     m0   : Total wave energy density (=ETOT)
!     m2   : Second moment of the variance spectrum (=ETOT2)
!     m4   : Fourth moment of the variance spectrum (=ETOT4)
!     k    : Wave number (=KWAVE(IS,1))
!     k~   : Mean wave number
!     Qb   : Fraction of breaking waves
!     sig  : Frequency (=SPCSIG(IS))
!     sig0 : Average zero crossing frequency
!     sig~ : Mean frequency
!     S~   : Overall steepness (STP_OV)
!     S_PM : Overall steepness for a Pierson-Moskowitz spectrum
!     th   : direction theta (=SPCDIR(ID))
!
!  4. Argument variables
!
!     AC2   : Action density
!     ACTOT : Total action density per gridpoint
!     CGO   : Group velocity (excluding current!)
!     DEP2  : Array containing water-depth
!     DISSC0: Dissipation coefficient as explicit part (meant for output)
!     DISSC1: Dissipation coefficient as implicit part (meant for output)
!     ETOT  : Total wave energy density
!     IDCMIN: Counter that indicates the minimum direction that is propagated in the sweep
!     IDCMAX: Counter that indicates the maximum direction that is propagated in the sweep
!     IMATDA: The values at the diagonal of the matrix that is solved numerically
!     IMATRA: The values at the right-hand side of the equation that is solved numerically
!     ISSTOP: Maximum counter in frequency space that is propagated within a sweep
!     KWAVE : Wavenumber
!     PLWCAP: Array containing the whitecapping source term for test-output
!     SPCDIR: (*,1); spectral directions (radians)
!             (*,2); cosine of spectral directions
!             (*,3); sine of spectral directions
!             (*,4); cosine^2 of spectral directions
!             (*,5); cosine*sine of spectral directions
!             (*,6); sine^2 of spectral directions
!     SPCSIG: Relative frequencies in computational domain in sigma-space
!     UFRIC : wind friction velocity
!
      INTEGER, INTENT(IN) :: ISSTOP, IDCMIN(MSC), IDCMAX(MSC)
!
      REAL, INTENT(IN)    :: AC2(MDC,MSC,MCGRD), DEP2(MCGRD)
      REAL, INTENT(IN)    :: ETOT
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL, INTENT(IN)    :: KWAVE(MSC,MICMAX)                            40.22
      REAL, INTENT(IN)    :: SPCDIR(MDC,6), SPCSIG(MSC)
      REAL, INTENT(OUT)   :: PLWCAP(MDC,MSC,NPTST)
      REAL, INTENT(INOUT) :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
      REAL, INTENT(INOUT) :: DISSC0(MDC,MSC,MDISP),DISSC1(MDC,MSC,MDISP)  40.67 40.12
      REAL, INTENT(IN)    :: UFRIC                                        40.53
      REAL, INTENT(IN)    :: CGO(MSC,MICMAX)                              40.53
!
!  6. Local variables
!
!     A     : Exponential term in the Longuet Higgins expression
!     B     : Value of saturation spectrum at current wavenumber
!     C_BJ  : Whitecapping coefficient according to Battjes-Janssen
!     C_K   : Whitecapping coefficient according to Komen
!     C_LH  : Whitecapping coefficient according to Longuet Higgins
!     EF    : Energy density spectrum in frequency domain
!             (azimuth-integrated frequency spectrum)
!     HM    : Maximum waveheight as used in the Battjes-Janssen expression
!     HRMS  : Significant wave height, based on total energy
!     ID    : Counter in directional space
!     ID1   : Counter in directional space
!     ID2   : Another counter in directional space
!     IDDUM : Counter in directional space within sweep limits
!     IENT  : Number of entries in this subroutine
!     IS    : Counter in frequency space
!     N1    : Exponent for the wavenumber term in the Komen expression
!     N2    : Exponent for the steepness term in the Komen expression
!     P     : Exponent of the relative saturation (B/Br)
!     QB_WC : The fraction of whitecapping waves in the Battjes-Janssen expression
!     SIG0  : Average zero-crossing frequency used in the Longuet Higgins expression
!     STP_OV: Overall steepness
!     STP_PM: Overall steepness for a Pierson-Moskowitz spectrum
!     WCAP  : Whitecapping source-term
!     WCIMPL: Implicit part of the whitecapping source-term
!
      INTEGER, SAVE     :: IENT = 0
      INTEGER           :: ID, IDDUM, IS, ID1, ID2
!
      REAL              :: A, C_BJ, HM, HRMS, N1, N2
      REAL              :: QB_WC, SIG0, STP_OV, STP_PM
!
      REAL, ALLOCATABLE :: C_K(:), C_LH(:), WCAP(:), WCIMPL(:)

      REAL              :: B, P, BRKD, FBR                                40.63 40.53
      REAL              :: PWCAP1, PWCAP2, PWCAP9, PWCAP10, PWCAP11       40.63
      REAL              :: EF(MSC)                                        40.53
!
!     8. REMARKS
!
!        alpha : PWCAP( 7)
!        C1    : PWCAP( 1)
!        C3    : PWCAP( 5)
!        delta : PWCAP(10)
!
!     9. STRUCTURE
!
!     ------------------------------------------------------------
!     Initialisation
!     ------------------------------------------------------------
!     Calculate needed parameters
!     ------------------------------------------------------------
!     If IWCAP = 1, 2, or 5; Calculate C_K
!     If IWCAP = 4 or 5; Calcualte C_BJ
!     If IWCAP = 3; Calculate C_LH
!     IF IWCAP = 7; Calculate with Alves and Banner (2003)
!     ------------------------------------------------------------
!     For frequency dependent part of the spectrum
!       Calculate dissipation term due to whitecapping
!     ------------------------------------------------------------
!     For the whole frequency domain
!       Fill the matrices and PLWCAP-array
!     ------------------------------------------------------------
!     End of SWCAP
!     ------------------------------------------------------------
!
! 13. Source text
!
      IF (LTRACE) CALL STRACE (IENT,'WCAP')
!
! Initialisation
!
      IF (ETOT.LE.0.) RETURN
      IF (ETOT2.LE.0.) RETURN
      IF (ETOT4.LE.0.) RETURN
      IF (ACTOT.LE.0.) RETURN
      IF (EDRKTOT.LE.0.) RETURN
!
      ALLOCATE (C_K(MSC), C_LH(MSC), WCAP(MSC), WCIMPL(MSC))
      WCIMPL(1:MSC) = 0.
!
! Calculate coefficients
!
      IF ((IWCAP.EQ.1).OR.
     &    (IWCAP.EQ.2).OR.
     &    (IWCAP.EQ.5)    ) THEN                                          40.30
!
! Calculate C_K
!
        STP_OV = KM_WAM * SQRT(ETOT)
        STP_PM = SQRT(PWCAP(2))
        N1     = PWCAP(11)
        N2     = 2. * PWCAP(9)
        C_K(:) = PWCAP(1) * (1. - PWCAP(10) +
     &           PWCAP(10) * (KWAVE(:,1) / KM_WAM)**N1) *
     &           (STP_OV / STP_PM)**N2
!
      ENDIF
!
      IF ((IWCAP.EQ.4).OR.
     &    (IWCAP.EQ.5)    ) THEN
!
! Calculate values for Hm and Qb
!
        HRMS   = SQRT(8. * ETOT)
        IF (IWCAP.EQ.4) HM = PWCAP(6) / KM01
        IF (IWCAP.EQ.5) HM = PWCAP(6) / (PWCAP(8) * KM_WAM)
        CALL FRABRE(HM, ETOT, QB_WC)
!
! Calculate C_BJ
!
        IF (HRMS.GE.HM) THEN
          C_BJ = PWCAP(7)  /  PI
        ELSE IF (HRMS.GT.0.) THEN
          C_BJ = (PWCAP(7) *  HM**2 * QB_WC) / (PI * HRMS**2)
        ELSE
          C_BJ = 0.
        END IF
      ENDIF
!
      IF (IWCAP.EQ.3) THEN
!
! Calculate C_LH
!
        SIG0 = SQRT(ETOT2 / ETOT)
!
!       A = -(1./8.)*(ETOT2**2/(ETOT*ETOT4))*(GRAV**2/(ETOT*SIG0**4))
!       rewrite to prevent underflow
!
        A = -(1./8.) * GRAV**2 / ETOT4
        DO IS=1, ISSTOP
!          C_LH(IS) = PWCAP(5) * SQRT((ETOT * SIG0**4) / GRAV**2) *
!     &               EXP(A) * SIG0 * (SPCSIG(IS) / SIG0)**2
!          rewrite to prevent underflow:
!
          C_LH(IS) = PWCAP(5) * EXP(A) * SQRT(ETOT2) * SPCSIG(IS)**2 /
     &               GRAV
        END DO
      END IF
!
! Calculate dissipation according to Alves & Banner (2003)                40.53
!
      IF ( IWCAP.EQ.7 ) THEN                                              40.53
!
! Calculate C_K                                                           40.63
!
! Note: use the default parameters of Komen et al. (1984) except Cds
!       which is slightly larger
!
        PWCAP1  = 3.00E-5                                                 40.63
        PWCAP2  = 3.02E-3                                                 40.63
        PWCAP9  = 2.                                                      40.63
        PWCAP10 = 0.                                                      40.63
        PWCAP11 = 1.                                                      40.63
!
        C_K    = 0.                                                       40.63
        STP_OV = KM_WAM * SQRT(ETOT)                                      40.63
        STP_PM = SQRT(PWCAP2)                                             40.63
        N1     = PWCAP11                                                  40.63
        N2     = 2. * PWCAP9                                              40.63
        C_K(:) = PWCAP1 * (1. - PWCAP10 +                                 40.63
     &           PWCAP10 * (KWAVE(:,1) / KM_WAM)**N1) *                   40.63
     &           (STP_OV / STP_PM)**N2                                    40.63
!                                                                         40.53
!  Loop to calculate B(k)                                                 40.53
!                                                                         40.53
        DO IS = 1, ISSTOP                                                 40.53
!                                                                         40.53
!  Calculate E(f)                                                         40.53
!                                                                         40.53
          EF(IS) = 0.                                                     40.53
          DO ID = 1,MDC                                                   40.53
            EF(IS) = EF(IS) + AC2(ID,IS,KCGRD(1))*SPCSIG(IS)*PI2*DDIR     40.53
          END DO                                                          40.53
!                                                                         40.53
!  Calculate saturation spectrum B(k) from E(f)                           40.53
!                                                                         40.53
          B = (1./PI2) * CGO(IS,1) * KWAVE(IS,1)**3 * EF(IS)              40.53
!                                                                         40.53
!  Calculate exponent P of the relative saturation (B/Br)                 40.53
!                                                                         40.53
          PWCAP(10)= 3. + TANH(25.76*(UFRIC*KWAVE(IS,1)/SPCSIG(IS)-0.1))  40.53
          BRKD = PWCAP(12)*(GRAV*KWAVE(IS,1)/(SPCSIG(IS)**2))**0          40.63
          P = PWCAP(10)                                                   40.63 40.53
          FBR = 0.5*(1. + TANH( 10.*( (B/BRKD)**0.5 - 1.)))               40.63
!                                                                         40.53
!  Calculate WCAP(IS) from B(k) and P                                     40.53
!                                                                         40.53
          STP_OV = KM_WAM * SQRT(ETOT)                                    40.53
          WCAP(IS) = PWCAP(1)*FBR*(B/BRKD)**(P/2.) *                      40.63 40.53
     &    STP_OV**PWCAP(9) * (KWAVE(IS,1)/KM_WAM)**PWCAP(11) *            40.53
     &    (GRAV**(0.5)*KWAVE(IS,1)**(0.5)/SPCSIG(IS))**(PWCAP(10)/2-1) *  40.53
     &    GRAV**(0.5)*KWAVE(IS,1)**(0.5)                                  40.53
     &    + (1.-FBR) * C_K(IS) * SIGM_10 * (KWAVE(IS,1) / KM_WAM)         40.63
                                                                          40.53
        END DO                                                            40.53
      END IF                                                              40.53
!
      IF ( IWCAP.LT.7 ) THEN                                              40.51 40.30
!
! Calculate the whitecapping source term WCAP(IS)
!
        DO IS=1, ISSTOP
          IF ((IWCAP.EQ.1).OR.
     &        (IWCAP.EQ.2).OR.
     &       ((IWCAP.EQ.5).AND.(C_BJ.LE.C_K(IS)))) THEN
            WCAP(IS) = C_K(IS) * SIGM_10 * (KWAVE(IS,1) / KM_WAM)
          ELSE IF (IWCAP.EQ.3) THEN
            WCAP(IS) = C_LH(IS)
          ELSE IF ((IWCAP.EQ.4).OR.
     &       ((IWCAP.EQ.5).AND.(C_BJ.GE.C_K(IS)))) THEN
            IF (IWCAP.EQ.4) WCAP(IS) = C_BJ*SIGM01 *(KWAVE(IS,1)/KM01  )
            IF (IWCAP.EQ.5) WCAP(IS) = C_BJ*SIGM_10*(KWAVE(IS,1)/KM_WAM)
!
! Calculate a term that is added to both sides of the equation to compensate
! for the strong non-linearity in the fraction of breaking waves Qb
!
            IF (HRMS.LT.HM) THEN
              WCIMPL(IS)=WCAP(IS) * ((1.-QB_WC)/((HRMS**2/HM**2)-QB_WC))
              WCAP(IS)  =WCAP(IS) + WCIMPL(IS)
            END IF
          ELSE
            CALL MSGERR(2,'Whitecapping is inactive')
            WRITE (PRINTF,*) 'Occurs in gridpoint: ', KCGRD(1)
          END IF
        END DO

      END IF
!
! Fill the diagonal of the matrix and the PLWCAP-array
!
      DO IS=1, ISSTOP
!
!        Only fill the values for the current sweep
!
         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD(IDDUM - 1 + MDC, MDC) + 1
            IMATDA(ID,IS)   = IMATDA(ID,IS)   + WCAP(IS)
            DISSC1(ID,IS,1) = DISSC1(ID,IS,1) + WCAP(IS)                  40.67 40.12
            IF (TESTFL) PLWCAP(ID,IS,IPTST) = -1.*(WCAP(IS)-WCIMPL(IS))
         END DO
      END DO
!
! Add the implicit part to the right-hand side, if appropriate
!
      IF ((IWCAP.EQ.4).OR.
     &    (IWCAP.EQ.5)) THEN
        DO IS=1, ISSTOP
!
!       Only fill the values for the current sweep
!
          DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD(IDDUM - 1 + MDC, MDC) + 1
            IMATRA(ID,IS)   = IMATRA(ID,IS) +
     &                        WCIMPL(IS) * AC2(ID,IS,KCGRD(1))
            DISSC0(ID,IS,1) = DISSC0(ID,IS,1) +                           40.67 40.12
     &                        WCIMPL(IS) * AC2(ID,IS,KCGRD(1))            40.12
          END DO
        END DO
      END IF
!
      DEALLOCATE (C_K, C_LH, WCAP, WCIMPL)
!
      RETURN
      END SUBROUTINE SWCAP
!****************************************************************
!
      SUBROUTINE BRKPAR (BRCOEF  ,ECOS    ,ESIN    ,AC2     ,             40.22
     &                   SPCSIG  ,DEP2    ,RDX     ,RDY     ,             41.03 30.72
     &                   KWAVE                              )             41.03
!
!****************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE                                                       40.22
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
!     40.02: IJsbrand Haagsma
!     40.22: Nico Booij
!     40.41: Marcel Zijlema
!     41.03: Andre van der Westhuysen
!
!  1. Updates
!
!            Jan. 97: New subroutine (Roeland Ris)
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.02, Oct. 00: KWAVE removed
!     40.22, Oct. 01: PSURF(2) is kept constant, BRCOEF added as argument
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.03, Feb. 09: extension with spatially varying breaker parameter
!                     according to Ruessink et al (2003)
!
!  2. Purpose
!
!     Determine the bottom slope in upwave direction and calculate
!     the slope dependent breaking parameter according to Nelson (1987)
!     Note that Nelson (1987) is used here since in Nelson (1994a,1994b)
!     an error is present in the equation.
!
!  3. Method
!
!     The breaker parameter is given by:
!
!     Hm / d =  0.55 + 0.88 exp ( -0.012 * cot beta)
!
!     with beta the angle the bed makes with the horizontal. This
!     above equation is only valid for positive slopes (Negative
!     slopes were not considered by Nelson. For very steep slopes
!     (>0.05 say) a very large breaker parameter is obtained (>>1).
!
!     To ensure wave breaking in laboratory cases (with very steep
!     slopes an upper limit of 0.81 (which corresponds to a bottom
!     slope of 0.01) is imposed on the model of Nelson.
!
!     For negative bottom slopes (not considered by Nelson) a value
!     of 0.73 is imposed (which is the average value in Table 2 of
!     Battjes and Janssen (1978).
!
!
!     Alternatively, the breaker index can be computed according to
!     Ruessink et al. (2003):
!
!     Hm / d = 0.76 kp * d + 0.29
!
!     with kp the peak wave number
!
!
!  4. Argument variables
!
      REAL, INTENT(OUT) :: BRCOEF    ! variable breaker coefficient       40.22

!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL, INTENT(IN)  :: SPCSIG(MSC)                                    30.72
      REAL, INTENT(IN)  :: AC2(MDC,MSC,MCGRD)  ! action densities         40.22
      REAL, INTENT(IN)  :: ECOS(MDC), ESIN(MDC)  ! Cos and Sin of Theta   40.22
      REAL, INTENT(IN)  :: DEP2(MCGRD)         ! depths at grid points    40.22
!     RDX, RDY:  coefficients to obtain spatial derivatives               40.22
      REAL, INTENT(IN)  :: RDX(10), RDY(10)                               40.08
      REAL, INTENT(IN)  :: KWAVE(MSC,MICMAX)                              41.03
!
!        INTEGERS :
!        ----------
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MCGRD       Maximum counter in geographical space
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!        KCGRD       Grid counter in central gridpoint (geo. space)
!
!
!        REALS:
!        ------
!        ETOTS       Total wave energy density in a particular
!                    direction (energy in tail neglected)
!        BRKVAR      Breaking coefficient
!        DIRDEG      Mean propagation direction of wave energy in
!                    degrees
!        DIRRAD      Mean propagation direction of wave energy in
!                    radians
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2         Action density
!        ECOS/ESIN   Cos. , sin of angle
!        DEP2        Depth
!        PSURF       Coefficients for breaking module
!
!  7. Common blocks used
!
!
!     5. SUBROUTINES CALLING
!
!        SINTGRL
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
!        ---
!
!     9. STRUCTURE
!
!     ------------------------------------------------------------
!     Calculate total energy per direction for all frequencies
!       determine action density per direction weighted with cos/sin
!       determine mean propagation direction of energy
!     ------------------------------------------------------------
!     calculate the depth derivative in the mean wave direction
!      according to dd/ds (see also subroutine SPROSD)
!     calculate the slope dependend breaking coefficient
!     ------------------------------------------------------------
!     End of NELSON
!     -------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER :: ID    ,IS      ! counters                                40.22
      INTEGER :: ISIGM                                                    41.03
!
!
      REAL  :: ETOTS ,EEX   ,EEY   ,
     &         EAD   ,SIGMA1,COSDIR,SINDIR,DDDX  ,                        40.22
     &         DDDY  ,DDDS  ,DETOT                                        40.22
      REAL  :: EMAX, ETD, KP, KPD                                         41.03
!
      INTEGER, SAVE :: IENT=0
      IF (LTRACE) CALL STRACE (IENT,'BRKPAR')
!
      IF ( ISURF.EQ.2 ) THEN                                              41.03
!
!        *** determine the average wave direction ***
!
         EEX   = 0.
         EEY   = 0.
         ETOTS = 0.
         DO ID = 1, MDC
            EAD = 0.
            DO IS = 1, MSC
               SIGMA1 = SPCSIG(IS)                                        30.72
               DETOT  = SIGMA1**2 * AC2(ID,IS,KCGRD(1))
               EAD    = EAD + DETOT
            ENDDO
            ETOTS = ETOTS + EAD
            EEX   = EEX + EAD * ECOS(ID)
            EEY   = EEY + EAD * ESIN(ID)
         ENDDO
!
         IF ( ETOTS .GT. 0.) THEN
            COSDIR = EEX / ETOTS
            SINDIR = EEY / ETOTS
         ELSE
            COSDIR = 1.
            SINDIR = 0.
         ENDIF
!
!        *** Determine bottom slope in average wave propagation direction ***
!
         DDDX =  RDX(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
     &         + RDX(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
         DDDY =  RDY(1) * (DEP2(KCGRD(1)) - DEP2(KCGRD(2)))
     &         + RDY(2) * (DEP2(KCGRD(1)) - DEP2(KCGRD(3)))
!
         DDDS = -1. * ( DDDX * COSDIR + DDDY * SINDIR )
!
!        *** calculate breaking coefficient according to Nelson (1987) ***
!
         IF ( DDDS .GE. 0. ) THEN
            DDDS   = MAX ( 1.E-6 , DDDS)
            BRCOEF = PSURF(4) + PSURF(7) * EXP ( -PSURF(8) / DDDS )       40.22
            BRCOEF = MIN ( PSURF(5) , BRCOEF )                            40.22
         ELSE
            BRCOEF = PSURF(6)                                             40.22
         ENDIF
!
      ELSE IF ( ISURF.EQ.3 ) THEN                                         41.03
!
!       calculate breaker index according to Ruessink et al (2003)
!
        EMAX = 0.
        ISIGM = -1
        DO IS = 1, MSC
           ETD = 0.
           DO ID = 1, MDC
              ETD = ETD + SPCSIG(IS)*AC2(ID,IS,KCGRD(1))*DDIR
           ENDDO
           IF (ETD.GT.EMAX) THEN
              EMAX  = ETD
              ISIGM = IS
           ENDIF
        ENDDO
        IF (ISIGM.GT.0) THEN
           KP = KWAVE(ISIGM,1)
        ELSE
           KP = 0.
        ENDIF
!
        KPD = KP*DEP2(KCGRD(1))
!
        IF ( KPD.LT.0.) THEN
           BRCOEF = 0.73
        ELSE
           BRCOEF = PSURF(4)*KPD + PSURF(5)
           BRCOEF = MIN( 1.2, BRCOEF)
           BRCOEF = MAX( 0.3, BRCOEF)
        ENDIF
!
      ENDIF
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 40 ) THEN
        WRITE(PRINTF,600) KCGRD(1), ATAN2(SINDIR,COSDIR)*180./PI,
     &                    DEP2(KCGRD(1)), DDDS, BRCOEF                    40.22
 600    FORMAT (' BRKPAR: point nr, dir, depth, slope, br.coeff:',
     &          I4,4(1X,E12.4))
      END IF
!
      RETURN
      END subroutine BRKPAR
!
!********************************************************************
!
      SUBROUTINE PLTSRC (PLWNDA        ,PLWNDB        ,
     &                   PLWCAP        ,PLBTFR        ,
     &                   PLWBRK        ,PLNL4S        ,
     &                   PLNL4D        ,PLTRI         ,
     &                   PLVEGT        ,                                  40.55
     &                   AC2           ,SPCSIG        ,                   40.00
     &                   DEP2          ,XYTST         ,
     &                                  KGRPNT        )                   40.00
!
!****************************************************************
!
      USE SWCOMM1                                                         40.41
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
!     1. UPDATE
!
!        40.00, Sep. 98: subroutine modified for new spectral file def.
!        40.00, Apr. 99: factor 2*PI added in 1d spectra
!        40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!     2. PURPOSE
!
!        store the source terms for the TESTFL gridpoint in a file
!
!     3. METHOD
!
!        ----
!
!     4. PARAMETERLIST
!
!        INTEGER
!        -------
!        MXC,MYC  Maximum counters in geographical space
!
!        REAL
!        ----
!
!        ARRAYS
!        -------
!        PLWNDA   2D    SWINEA
!        PLWNDB   2D    SWINEB
!        PLWCAP   2D    SWCAPE
!        PLBTFR   2D    SWBOTE
!        PLVEGT   2D    SVEGET
!        PLWBRK   2D    SWSURF
!        PLNL4    2D    SWNL
!        PLTRI    2D    TRIADS
!
!     5. subroutines calling
!
!
!     6. subroutines used
!
!        WRSPEC
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ----------------------------------------------------------------
!   If Mode Timedep is active
!   Then write time
!   Else write iteration
!   ----------------------------------------------------------------
!   For all test points do
!       If IFS1D > 0
!       Then For all spectral frequencies do
!                integrate action and source terms over direction
!                write action and source terms to file IFS1D
!       ------------------------------------------------------------
!       If IFS2D > 0
!       Then write action densities and source terms to file IFS2D
!       ------------------------------------------------------------
!       Set source terms equal to 0
!   ----------------------------------------------------------------
!
!     10. SOURCE
!
      INTEGER     IS    ,ID
!
!     SIGACT      product of sigma and action density, i.e. energy density
!     WCAP        integral of whitecapping dissipation
!
      REAL        WCAP  ,BTFR  ,WBRK  ,NL4   ,FAC   ,SIGACT,
     &            VEGT  ,                                                 40.55
     &            NL4S  ,NL4D  ,TRIA  ,ENERGY,ENRSIG
!
      INTEGER     XYTST(*),KGRPNT(MXC,MYC)                                40.80 30.21
!
!
      REAL        AC2(MDC,MSC,MCGRD)          ,
     &            SPCSIG(MSC)                 ,                           40.00
     &            PLWNDA(MDC,MSC,NPTST)       ,
     &            PLWNDB(MDC,MSC,NPTST)       ,
     &            PLWCAP(MDC,MSC,NPTST)       ,
     &            PLBTFR(MDC,MSC,NPTST)       ,
     &            PLVEGT(MDC,MSC,NPTST)       ,                           40.55
     &            PLWBRK(MDC,MSC,NPTST)       ,
     &            PLNL4S(MDC,MSC,NPTST)       ,
     &            PLNL4D(MDC,MSC,NPTST)       ,
     &            PLTRI (MDC,MSC,NPTST)       ,
     &            DEP2(MCGRD)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'PLTSRC')
!
!     *** compute the 1D spectra ***
!
      DO 300 IPTST = 1, NPTST
        IF (OPTG.NE.5) THEN                                               40.80
           LXDMP = XYTST(2*IPTST-1)
           LYDMP = XYTST(2*IPTST)
           INDX  = KGRPNT(LXDMP,LYDMP)
        ELSE                                                              40.80
           INDX = XYTST(IPTST)                                            40.80
        ENDIF                                                             40.80
        IF (IFPAR.GT.0) THEN
!
!         *** for the parameter output we first integrate ***
!         *** over all frequencies and directions         ***             40.00
!
          ENERGY = 0.
          ENRSIG = 0.
          SWND = 0.
          WCAP = 0.
          BTFR = 0.
          VEGT = 0.
          WBRK = 0.
          NL4  = 0.
          TRIA = 0.
          DO 60 IS = 1, MSC
            DO 50 ID = 1, MDC
!
!             *** ENERGY density ***
!
              SIG2AC = SPCSIG(IS)**2 * AC2(ID,IS,INDX)
              ENERGY = ENERGY + SIG2AC                                    40.00
              ENRSIG = ENRSIG + SPCSIG(IS) * SIG2AC                       40.00
!
!             *** wind input ***
!
              SWND = SWND + SPCSIG(IS)**2 * PLWNDA(ID,IS,IPTST) +
     &                 PLWNDB(ID,IS,IPTST) * SIG2AC                       40.00
!
!             *** dissipation processes ***
!
              WCAP = WCAP + PLWCAP(ID,IS,IPTST) * SIG2AC                  40.00
              BTFR = BTFR + PLBTFR(ID,IS,IPTST) * SIG2AC
              VEGT = VEGT + PLVEGT(ID,IS,IPTST) * SIG2AC                  40.55
              WBRK = WBRK + PLWBRK(ID,IS,IPTST) * SIG2AC
!
!             *** nonlinear interactions ***
!
              TRIA = TRIA + ABS(PLTRI (ID,IS,IPTST)) * SPCSIG(IS)**2      40.85 40.00
!
              IF ( IQUAD .EQ. 1) THEN
                NL4  = NL4  + ABS(PLNL4D(ID,IS,IPTST) * SIG2AC +
     &                            PLNL4S(ID,IS,IPTST) * SPCSIG(IS)**2)    40.85 40.00
              ELSE
                NL4  = NL4  + ABS(PLNL4S(ID,IS,IPTST)) * SPCSIG(IS)**2    40.85 40.00
              END IF
!
  50        CONTINUE
  60      CONTINUE
!
          IF (ENERGY.GT.0.) THEN                                          40.00
            ENERGY = ENERGY * FRINTF * DDIR
            ENRSIG = ENRSIG * FRINTF * DDIR
            SWND   = SWND   * FRINTF * DDIR
            WCAP   = WCAP   * FRINTF * DDIR
            BTFR   = BTFR   * FRINTF * DDIR
            VEGT   = VEGT   * FRINTF * DDIR                               40.55
            WBRK   = WBRK   * FRINTF * DDIR
            TRIA   = TRIA   * FRINTF * DDIR
            NL4    = NL4    * FRINTF * DDIR
!
            WRITE (IFPAR, 70) 4.*SQRT(ENERGY), PI2*ENERGY/ENRSIG,
     &                          SWND, WCAP, BTFR, VEGT, WBRK, TRIA, NL4   40.55 40.00
  70        FORMAT(9(1X,E12.4))
          ELSE
            WRITE (IFPAR, 70) OVEXCV(10), OVEXCV(28), OVEXCV(7),          40.41
     &          OVEXCV(7),                                                40.55
     &          OVEXCV(7), OVEXCV(7), OVEXCV(7), OVEXCV(7), OVEXCV(7)     40.00
          ENDIF
        ENDIF
!
        IF (IFS1D.GT.0) THEN
          IF (DEP2(INDX).LE.0.) THEN
            WRITE (IFS1D, 80) 'NODATA'                                    40.00
          ELSE
            WRITE (IFS1D, 80) 'Test point ', IPTST                        40.00
  80        FORMAT (A, I6)                                                40.00
!
!         *** for the output of 1D spectra we integrate over ***
!         *** all directions                                 ***
!
            DO 160 IS = 1, MSC
              SWND = 0.
              WCAP = 0.
              BTFR = 0.
              VEGT = 0.
              WBRK = 0.
              NL4S = 0.
              NL4D = 0.
              TRIA = 0.
              ENERGY = 0.
              DO 150 ID = 1, MDC
                SIGACT = SPCSIG(IS) * AC2(ID,IS,INDX)
!
!               *** wind input ***
!
                SWND = SWND + PLWNDA(ID,IS,IPTST) * SPCSIG(IS) +
     &                        PLWNDB(ID,IS,IPTST) * SIGACT                40.00
!
!               *** dissipation processes ***
!
                WCAP = WCAP + PLWCAP(ID,IS,IPTST) * SIGACT                40.00
                BTFR = BTFR + PLBTFR(ID,IS,IPTST) * SIGACT                40.00
                VEGT = VEGT + PLVEGT(ID,IS,IPTST) * SIGACT                40.55
                WBRK = WBRK + PLWBRK(ID,IS,IPTST) * SIGACT                40.00
!
!               *** nonlinear interactions ***
!
                TRIA = TRIA + PLTRI (ID,IS,IPTST) * SPCSIG(IS)            40.00
!
                NL4S = NL4S + PLNL4S(ID,IS,IPTST) * SPCSIG(IS)            40.00
!
                IF ( IQUAD .EQ. 1) THEN
                  NL4D = NL4D + PLNL4D(ID,IS,IPTST) * SIGACT              40.00
                END IF
!
!               *** energy density ***
!
                ENERGY = ENERGY + SIGACT                                  40.00

 150          CONTINUE
              NL4 = NL4S + NL4D
!             factor 2*PI introduced to account for density per Hz        40.00
!             instead of per rad/s.
!             factor DDIR is due to integration over directions
              FAC = PI2 * DDIR
              WRITE (IFS1D,170) ENERGY*FAC, SWND*FAC, WCAP*FAC,           40.00
     &                          BTFR*FAC, VEGT*FAC, WBRK*FAC, TRIA*FAC,   40.55 40.00
     &                          NL4*FAC                                   40.00
 170          FORMAT(9(1X,E12.4))                                         40.00
 160        CONTINUE
          ENDIF
        ENDIF
!
!       output of 2D distributions of source terms
!
        IF (IFS2D.GT.0) THEN
          IF (DEP2(INDX).LE.0.) THEN
            DO LOOP = 1, 7
              WRITE (IFS2D, 80) 'NODATA'                                  40.03
            ENDDO
          ELSE
            DO IS = 1, MSC
              DO ID = 1, MDC
                SIGACT = SPCSIG(IS) * AC2(ID,IS,INDX)
                PLWNDA(ID,IS,IPTST) = PLWNDA(ID,IS,IPTST) * SPCSIG(IS)    40.00
     &                              + PLWNDB(ID,IS,IPTST) * SIGACT
                PLWCAP(ID,IS,IPTST) = PLWCAP(ID,IS,IPTST) * SIGACT
                PLBTFR(ID,IS,IPTST) = PLBTFR(ID,IS,IPTST) * SIGACT
                PLVEGT(ID,IS,IPTST) = PLVEGT(ID,IS,IPTST) * SIGACT        40.55
                PLWBRK(ID,IS,IPTST) = PLWBRK(ID,IS,IPTST) * SIGACT
                PLNL4S(ID,IS,IPTST) = PLNL4S(ID,IS,IPTST) * SPCSIG(IS)
     &                              + PLNL4D(ID,IS,IPTST) * SIGACT
!               PLWNDB is used temporarily for energy density             40.00
                PLWNDB(ID,IS,IPTST) = SIGACT
              ENDDO
            ENDDO
            CALL WRSPEC (IFS2D, PLWNDB(1,1,IPTST))                        40.00
            CALL WRSPEC (IFS2D, PLWNDA(1,1,IPTST))
            CALL WRSPEC (IFS2D, PLWCAP(1,1,IPTST))
            CALL WRSPEC (IFS2D, PLBTFR(1,1,IPTST))
            CALL WRSPEC (IFS2D, PLVEGT(1,1,IPTST))                        40.55
            CALL WRSPEC (IFS2D, PLWBRK(1,1,IPTST))
            CALL WRSPEC (IFS2D, PLTRI (1,1,IPTST))
            CALL WRSPEC (IFS2D, PLNL4S(1,1,IPTST))
          ENDIF
        ENDIF
!
!       *** set arrays zero: ***
!
        DO 200 IS = 1, MSC
          DO 100 ID = 1, MDC
            PLWNDA(ID,IS,IPTST) = 0.
            PLWNDB(ID,IS,IPTST) = 0.
            PLWCAP(ID,IS,IPTST) = 0.
            PLBTFR(ID,IS,IPTST) = 0.
            PLVEGT(ID,IS,IPTST) = 0.
            PLWBRK(ID,IS,IPTST) = 0.
            PLTRI (ID,IS,IPTST) = 0.
            PLNL4S(ID,IS,IPTST) = 0.
            PLNL4D(ID,IS,IPTST) = 0.
 100      CONTINUE
 200    CONTINUE
!
 300  CONTINUE
!
      RETURN
!     end of subroutine PLTSRC
      END
