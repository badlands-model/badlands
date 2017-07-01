!
!     SWAN/COMPU    file 1 of 5
!
!
!     PROGRAM SWANCOM1.FOR
!
!     This subroutine SWANCOM1 of the main program SWAN
!     includes the next subroutines :
!
!     SWCOMP  (main subroutine for the computational module)              31.02
!     SWOMPU  (carries out computation for one grid point)                31.02
!     SWPRSET (print all the settings used in SWAN run)                   40.80
!     SACCUR  (calculate the accuracy and check if the iteration process
!              can be terminated)
!     INSAC   (initialize the values for the calculation of accuracy)     40.00
!     ACTION  (fill the arrays with the derivatives of the action eq.)
!     SINTGRL (calculate some general wave (integral) parameters)         40.02
!     SOLPRE  (preparation before solving the system)
!     SOLMAT  (solve tri-diagonal system in absence of a current)
!     SOLMT1  (solve tri-diagonal system in presence of a current;
!              almost identical to SOLMAT, however, space in theta can
!              be periodic)
!     SOURCE  (fill the array with the source terms)
!     PHILIM  (limit the change in action density between two iterations)
!     HJLIM   (limit the change in action density between two iterations  40.61
!              based on Hersbach and Janssen limiter)                     40.61
!     RESCALE (remove negative values from action density)
!     SWSIP   (solve penta-diagonal system in spectral space by means     40.23
!              of Stone's SIP solver)                                     40.23
!     SWSOR   (solve penta-diagonal system in spectral space by means     40.41
!              of point SOR method)                                       40.41
!     SWMTLB  (compute bounds of thread loop)                             40.31
!     SWSTPC  (calculate the accuracy and check if the iteration process  40.41
!              can be terminated based on curvature of Hs)                40.41
!     SETUPP  (compute the wave-induced setup for a one-dimensional and a 32.01
!              two-dimensional run. Note that the one-dimensional mode of 32.01
!              SWAN has been coded in this project (H3268))               32.01
!     SETUP2D (computation of the change of waterlevel by waves,          31.03
!              a 2D Poisson equation in general coordinates is solved)    31.03
!
!******************************************************************
!
      SUBROUTINE SWCOMP (AC1        ,AC2        ,                         40.22
     &                   COMPDA     ,                                     40.22
     &                   SPCDIR     ,SPCSIG     ,                         30.72
     &                   XYTST      ,
     &                   IT         ,KGRPNT     ,
     &                   XCGRID     ,YCGRID     ,                         30.72
     &                   CROSS      )                                     40.31 40.30
!
!******************************************************************

      USE TIMECOMM                                                        40.41
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_PARALL                                                        40.31
      USE m_constants                                                     40.41
      USE m_xnldata                                                       40.41
      USE m_fileio                                                        40.41
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
!     30.74: IJsbrand Haagsma (Include version)
!     30.75: IJsbrand Haagsma (Bug fix)
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     31.03: Annette Kieftenburg
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     33.08: W. Erick Rogers (some S&L scheme-related changes)
!     33.10: Nico Booij and Erick Rogers
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!     40.21: Agnieszka Herman
!     40.22: John Cazes and Tim Campbell
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.31: Tim Campbell and John Cazes
!     40.31: Andre van der Westhuysen
!     40.41: Andre van der Westhuysen
!     40.41: Marcel Zijlema
!     40.41: Andre van der Westhuysen
!
!  1. Updates
!
!     30.72, Nov. 97: Declaration of MSC4MI, MSC4MA, MDC4MI, MDC4MA and
!                     ISTAT removed because they are common and already
!                     declared in the INCLUDE file
!     30.72, Nov. 97: ITERMX can be chosen freely with NUM ACCUR also in dynamic
!                     mode. Default ITERMX=6. Needs extensive testing
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     32.02, Jan. 98: Introduced 1D-version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: call WINDP0 removed, function taken over by WINDP1
!     30.72, Mar. 98: Current switched off for the first iteration, when
!                     preconditining is required
!     30.72, Mar. 98: Writes the result of the iteration step to the PRINT
!                     file
!     30.75, Mar. 98: Renamed SLOW to SIGLOW, because SLOW was used only locally
!     31.03, Feb. 98: Call SETUPP added, initialisation of array SETPDA
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Jan. 99: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     33.08, July 98: some S&L scheme-related changes
!     30.82, July 99: Corrected argumentlist SETUPP
!     40.00, July 99: argument KQUAD removed from call PLTSRC
!     30.82, Sep. 99: Modified messages in case of non-convergence
!     33.10, Jan. 00: minor changes re: the SORDUP scheme
!     40.03, Mar. 00: Ursell number is now array (value for each grid point)
!     40.02, Oct. 00: Avoided real/int conflict by introducing replacing
!                     RWAREA for WAREA in FAC4WW and SETUPP
!     40.13, Mar. 01: comments changed;
!                     order of calling SWAPAR and SPROXY changed
!                     message concerning lack of convergence only to print file
!                     in nonstationary cases
!     40.21, Aug. 01: implementation of diffraction
!     40.22, Sep. 01: WAREA, LWAREA, and RWAREA structures removed        40.22
!                     and replaced with allocated arrays to ease          40.22
!                     OpenMP implementation.                              40.22
!     40.22, Sep. 01: OpenMP directives were added to parallelize the     40.22
!                     outer Y loop for the call to SWOMPU in the sweep    40.22
!                     across the computational grid.                      40.22
!     40.22, Sep. 01: Added logical array LLOCK for thread management     40.22
!                     during parallel operation.  It will not affect      40.22
!                     serial execution.                                   40.22
!     40.22, Sep. 01: Changed array definitions to use the parameter      40.22
!                     MICMAX instead of ICMAX.                            40.22
!     40.17, Dec. 01: Implemented Multiple DIA
!     40.23, Aug. 02: Print of CPU times added
!     40.23, Aug. 02: Print of use of limiter and rescaling
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Jul. 03: some improvements and corrections w.r.t. OpenMP
!     40.31, Sep. 03: Under-relaxation parameter is set to zero in first
!                     iteration because of the first guess
!     40.41, May  04: Implemented XNL (Webb-Resio-Tracy) method for
!                     quadruplet interactions
!     40.41, Jun. 04: Implementation of curvature-based convergence check
!     40.41, Aug. 04: some code optimization
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     The aim of this model is to simulate the wave energy in
!     shallow water areas. In the subroutine SWCOMP the main processes
!     taking place in the shallow water zone are determined in
!     several subroutines.
!     The input for this subroutine comes from SWANPRE1 and SWANPRE2.
!     The output is send to the subroutines SWANOUT1 and SWANOUT2.
!     The output consist of some characteristic
!     wave parameters and the wave action density. The equations are
!     all based on the action density N which is a function of the
!     spatial position (x,y), the relative frequency (s) and the
!     spectral direction (d).
!
!  3. Method
!
!     Keywords:
!     Action density, propagation terms, refraction, reflection,
!     white capping, wave breaking, bottom friction, nonlinear
!     and nonhomogeneous wind- and current-fields, wave blocking,
!     fully spectral description, nonlinear wave-wave interaction,
!     higher order upwind schemes, flux limiting, SIP solver
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
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
!
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!
!     INTEGERS:
!     --------------------------------------------------------------
!     IC                Dummy variable: ICode gridpoint:
!                       IC = 1  Top or Bottom gridpoint
!                       IC = 2  Left or Right gridpoint
!                       IC = 3  Central gridpoint
!                       Whether which value IC has, depends of the sweep
!                       If necessary, IC can be enlarged by increasing
!                       the array size of ICMAX
!     ITER              Counter of iterations per 4 sweeps for accuracy
!     ITERMX            Maximum number of iterations in model
!     IX                Counter of gridpoints in x-direction
!     IY                Counter of gridpoints in y-direction
!     IS                Counter of relative frequency band
!     IT                Counter in time space
!     ID                Counter of directional distribution
!     IBOT              Indicator for bottom friction
!                       IBOT = 0  no bottom friction dissipation
!                       IBOT = 1  Jonswap bottom dissipation model
!                       IBOT = 2  Dingemans bottom dissipation model
!                       IBOT = 3  Madsen bottom dissipation model
!     ICUR              Indicator for current
!     ISURF             Indicator for wave breaking
!     ITRIAD            Indicator for nonlinear triad interactions
!     IQUAD             Indicator for nonlinear quadruplet interactions
!     IWCAP             Indicator for white capping
!     IWIND             Indicator for which wind generation model is used
!                       IWIND = 1 first generation wind growth model
!                       IWIND = 2 second generation wind growth model
!                       IWIND = 3 third generation wind growth model
!     IREFR             indicator for refraction (can be tuned off)
!     ITFRE             indicator for transport of action in frequency
!                       space
!     ICMAX             Maximum array size for the points in the molecule
!     KSX      input    Dummy variable to get the right sign in the
!                       numerical difference scheme in X-direction
!                       depending on the sweep direction
!     KSY      input    Dummy variable to get the right sign in the
!                       numerical difference scheme in Y-direction
!                       depending on the sweep direction
!     MXC               Maximum counter of gridppoints in x-direction in
!                       computational model: (XLEN/DX + 1 )
!     MYC               Maximum counter of gridppoints in y-direction in
!                       computational model: (YLEN/DY + 1 )
!     MSC               Maximum counter of relative frequency in
!                       computational model
!     MDC               Maximum counter of directional distribution in
!                       computational model (2PI / DDIR + 1)
!     MTC               Maximum counter of the time, i.e.:
!                       (total time in proto type) / (time step)
!     MBOT              Maximum array size for PBOT
!     MSURF             Maximum array size for PSSURF
!     MTRIAD            Maximum array size for PTRIAD
!     MWCAP             Maximum array size for PWCAP
!     MWIND             Maximum array size for PWIND
!     MNISL             Minimum sigma-index occured in applying limiter
!     MXNFL             Maximum number of use of limiter in a spectral
!                       space
!     MXNFR             Maximum number of use of rescaling in a spectral
!                       space
!     NPFL              Number of geographical points in which limiter
!                       is used
!     NPFR              Number of geographical points in which rescaling
!                       is used
!     NWETP             Total number of wet gridpoints
!     IDEBUG            Level of debug output:
!                       0 = no output
!                       1 = print of statistics w.r.t. use of limiter
!                           and rescaling
!     I1GRD             Lower index for thread loop over spatial grid     40.31
!     I2GRD             Upper index for thread loop over spatial grid     40.31
!     I1MYC             Lower index for thread loop over y-grid row       40.31
!     I2MYC             Upper index for thread loop over y-grid row       40.31
!
!     REALS:
!     --------------------------------------------------------------
!
!     ALEN              Part of side length of an angle side
!     BETA              Angle between DX end DY
!     BLEN              Part of side length of an angle side
!     DIR               Spectral direction (i.e., ID*DDIR)
!     DX       input    Length of spatial cell in X-direction
!     DY       input    Length of spatial cell in Y-direction
!     DS       input    Width of frequency band (is not constant because
!                       of the logarithmic distribution of the frequency
!     DDIR     input    Width of directional band
!     DT       input    Time step
!     DDX      input    Same as DX but with correct sign depending of the
!                       direction of the sweep (+1. OR -1. ) no input
!     DDY      input    Same as DY but with correct sign depending of the
!                       direction of the sweep (+1. OR -1. ) no input
!     FAC_A             Factor representing the influence of the action-
!                       density depening of the propagation velocity
!     FAC_B             Factor representing the influence of the action-
!                       density depending of the propagation velocity
!     GAMMA             PI - alpha - beta
!     HM                Maximum waveheight (breaking source term)
!     GRAV     input    Gravitational acceleration
!     FRAC              Fraction of total wet points
!
!     one and more dimensional arrays:
!     ---------------------------------
!
!     AC1       4D    Action density as function of D,S,X,Y at time T
!     AC2       4D    (Nonstationary case) action density as function
!                     of D,S,X,Y at time T+DT
!     CGO       2D    Group velocity as function of IC and IS in the
!                     direction of wave propagation in absence of currents
!     CAX       3D    Wave transport velocity in X-direction, function of
!                     (ID,IS,IC)
!     CAY       3D    Wave transport velocity in Y-direction, function of
!                     (ID,IS,IC)
!     CAS       3D    Wave transport velocity in S-direction, function of
!                     (ID,IS,IC)
!     CAD       3D    Wave transport velocity in D-dirction, function of
!                     (ID,IS,IC)
!     COMPDA    3D    array containing depth and other arrays of (IX,IY)  20.39
!                     JDP1    Depth as function of X and Y at time T
!                     JDP2    (Nonstationary case) depth as function of X and Y
!                             at time T+DT
!                     JVX1    X-component of current velocity of X and Y
!                             at time T
!                     JVX2    (Nonstationary case) X-component of current
!                             velocity in (X,Y) at time T+DT
!                     JVY1    Y-component of current velocity in (X,Y)
!                             at time T
!                     JVY2    (Nonstationary case) Y-component of current
!                             velocity in (X,Y) at time T+DT
!                     JWX2    X-component of wind velocity in (X,Y)       40.00
!                             at time T+DT (nonstationary case)           40.00
!                     JWY2    Y-component of wind velocity in (X,Y)       40.00
!                             at time T+DT (nonstationary case)           40.00
!                     JUBOT   Absolute orbital velocity in a gridpoint (IX,IY)
!     SWTSDA    4D    intermediate data computed for the test points;
!                     there are MTSVAR subarrays:
!                     JPWNDA   wind source term part A
!                     JPWNDB   wind source term part B
!                     JPWCAP   whitecapping source term
!                     JPBTFR   bottom friction
!                     JPVEGT   vegetation dissipation
!                     JPWBRK   surf breaking
!                     JP4S     quadruplet interactions
!                     JP4D     quadruplet interactions
!                     JPTRI    triad interactions
!     ALIMW     1D    Maximum energy by wind growth. This dummy array is
!                     used because the maximum value has to be checked
!                     direct after the solver of the tri-diagonal matrix
!                     see the subroutine SOLMAT
!     GROWW     1D    Check for a certain frequency if the waves are
!                     growing or not in a spectral direction (LOGICAL)
!     HSAC0     2D    Represent the significant wave height at iter-2
!     HSAC1     2D    Represent the significant wave height at iter-1
!     HSAC2     2D    Represent the significant wave height at iter
!     HSDIFC    2D    Represent Hs(iter) - Hs(iter-2) meant for
!                     computation of curvature of Hs
!     IMATDA    2D    Coefficients of main diagonal of matrix
!     IMATLA    2D    Coefficients of lower diagonal of matrix in theta-space
!     IMATUA    2D    Coefficients of upper diagonal of matrix in theta-space
!     IMAT5L    2D    Coefficients of lower diagonal of matrix in sigma-space
!     IMAT6U    2D    Coefficients of upper diagonal of matrix in sigma-space
!     IMATRA    2D    Coefficients of right hand side
!     KWAVE     2D    wavenumber as function of the relative frequency S
!                     and position IC(ix,iy)
!     PBOT      1D    Coefficient for the bottom friction models
!     PSURF     1D    Coefficient for the wave breaking model
!     PTRIAD    1D    Coefficient for the triad interaction model
!     PWCAP     1D    Coefficient for the white capping model
!     PWIND     1D    Coefficient for the wind growth model
!     SACC0     2D    Represents the mean wave frequency at iter-2
!     SACC1     2D    Represents the mean wave frequency at iter-1
!     SACC2     2D    Represents the mean wave frequency at iter
!     TMDIFC    2D    Represent Tm(iter) - Tm(iter-2) meant for
!                     computation of curvature of Tm
!     PWTAIL    1D    coefficients for tail of spectrum
!     IDCMIN    1D    frequency dependent counter in directional space
!                     no current <---> current
!     IDCMAX    1D    frequency dependent counter in directional space
!                     no current <---> current
!     ISCMIN    1D    frequency dependent counter in frequency space
!                     no current <---> current
!     ISCMAX    1D    frequency dependent counter in frequency space
!                     no current <---> current
!     ANYBIN    2D    Set for a particular bin TRUE or FALSE depending on
!                     propagation velocities within a sweep
!     WWINT     1D    Counters for 4 wave-wave interactions
!     WWAWG     1D    Weight coefficients for the 4 wave-wave interactions
!     WWSWG     1D    Weights coefficients for the 4 wave-wave interactions
!                     for the semi-implicit computation
!     ISLMIN    1D    Lowest sigma-index occured in applying limiter
!     NFLIM     1D    Number of frequency use of limiter in each
!                     geographical point
!     NRSCAL    1D    Number of frequency use of rescaling in each
!                     geographical point
!     AC2LOC    2D    help array containing action density for
!                     sending/receiving with MPI
!     IARR      1D    help array of type integer for MPI communication
!     ARR       1D    help array of type real for MPI communication
!     REFLSO    2D    contribution to the source term due to reflection
!
!     Coefficients for the arrays:
!     -----------------------------
!                         default
!                         value:
!
!     PBOT(1)   = CFC      0.005    (Collins equation)
!     PBOT(2)   = CFW      0.01     (Collins equation)
!     PBOT(3)   = GAMJNS   0.067    (Jonswap equation)
!     PBOT(4)   = MF      -0.08     (Madsen equation)
!     PBOT(5)   = KN       0.05     (bottom roughness)
!
!     ISURF                1        (Constant breaking coefficient)
!                          2        (variable breaking coefficient
!                                    according to Nelson (1994))
!     PSURF(1)  = ALFA     1.0      (Battjes Janssen)
!     PSURF(2)  = GAMMA    0.73     (breaking criterium)
!
!     PWCAP(1)  = ALFAWC   2.36e-5  (Empirical coefficient)
!     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
!
!     PWIND(1)  = CF10     188.0    (second generation wind growth model)
!     PWIND(2)  = CF20     0.59     (second generation wind growth model)
!     PWIND(3)  = CF30     0.12     (second generation wind growth model)
!     PWIND(4)  = CF40     250.0    (second generation wind growth model)
!     PWIND(5)  = CF50     0.0023   (second generation wind growth model)
!     PWIND(6)  = CF60    -0.2233   (second generation wind growth model)
!     PWIND(7)  = CF70     0.       (second generation wind growth model)
!     PWIND(8)  = CF80    -0.56     (second generation wind growth model)
!     PWIND(9)  = RHOAW    0.00125  (density air / density water)
!     PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!     PWIND(11) = CDRAG    0.0012   (drag coefficient)
!     PWIND(12) = UMIN     1.0      (minimum wind velocity)
!     PWIND(13) = PMLM     0.13     (  )
!
!     PNUMS(1)  = DREL     relative error in Hs and Tm
!     PNUMS(2)  = DHABS    absolute error in Hs
!     PNUMS(3)  = DTABS    absolute error in Tm
!     PNUMS(4)  = NPNTS    number of points were accuracy is reached
!
!     PNUMS(5)  = NOT USED
!     PNUMS(6)  = CDD      blending parameter for finite differences
!                          in theta space
!     PNUMS(7)  = CSS      blending parameter for finite differences
!                          in sigma space
!     PNUMS(8)  = NUMFRE   numerical scheme in frequency space :
!                          1) implicit scheme
!                          2) explicit scheme CFL limited
!                          3) explicit scheme filter after iteration
!     PNUMS(9)  = DIFFC    if explicit scheme is used, then numerical
!                          diffusion coefficient can be chosen
!     PNUMS(12) = EPS2     termination criterion in relative sense for a
!                          penta-diagonal solver
!     PNUMS(13) = OUTP     request for output for a penta-diagonal solver
!     PNUMS(14) = NITER    maximum number of iterations for a penta-diagonal
!                          solver
!     PNUMS(15) = DHOVAL   global error in Hs
!               = CURVAT   curvature of Hs meant for convergence check
!     PNUMS(16) = DTOVAL   global error in Tm01
!     PNUMS(17) = CDLIM    coefficient of limitation of Ctheta
!     PNUMS(18) = FROUDMAX maximum Froude number for reduction of currents
!     PNUMS(19) = CFL      CFL criterion for option explicit scheme
!                          in frequency space (see PNUMS(8))
!     PNUMS(20) = GRWMX    maximum growth in spectral bin
!
!     PNUMS(21) = STOPC    type of stopping criterion:
!                          0: standard SWAN based on relative and global
!                             errors of Hs and Tm01,
!                          1: based on absolute, relative and curvature
!                             errors of Hs
!
!     PNUMS(30) = ALFA     relaxation parameter for under-relaxation method
!
!     arrays for the 4-wave interactions:
!
!     WWINT ( 1 = IDP    WWAWG ( = AGW1    WWSWG ( = SWG1
!             2 = IDP1           = AWG2            = SWG2
!             3 = IDM            = AWG3            = SWG3
!             4 = IDM1           = AWG4            = SWG4
!             5 = ISP            = AWG5            = SWG5
!             6 = ISP1           = AWG6            = SWG6
!             7 = ISM            = AWG7            = SWG7
!             8 = ISM1           = AWG8 )          = SWG8  )
!             9 = ISLOW
!             10= ISHGH
!             11= ISCLW
!             12= ISCHG
!             13= IDLOW
!             14= IDHGH
!             15= MSC4MI
!             16= MSC4MA
!             17= MDC4MI
!             18= MDC4MA
!             19= MSCMAX
!             20= MDCMAX
!             21= IDPP
!             22= IDMM
!             23= ISPP
!             24= ISMM )
!
!
!  6. Local variables
!
!     SIGLOW: recommended lowest frequency when TRIADS are activated
!
      REAL    SIGLOW
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     INSAC
!     SWOMPU
!     PLTSRC
!     SETUPP
!     SACCUR
!     SWSTPC
!TIMG!     SWTSTA                                                              40.23
!TIMG!     SWTSTO                                                              40.23
!     SWREDUCE                                                            40.30
!     SWRECVAC                                                            40.30
!     SWSENDAC                                                            40.30
!     MSGERR : Handles error messages according to severity               40.41
!     NUMSTR : Converts integer/real to string                            40.41
!     TXPBLA : Removes leading and trailing blanks in string              40.41
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWANPREn, SWANOUTn
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     SWCOMP is the main subroutine is and called of the main program SWAN.
!     The main program SWAN is build of three main subroutines:
!
!     1. SWREAD    (preparation of the computation (reading parameters))
!     2. SWCOMP    (computation of the action densities (discussed below))
!     3. SWOUTP    (output of the computation)
!
!     In this part the subroutine SWCOMP is discussed:
!
!     SWCOMP
!     ======
!      |
!      |------+  INSAC                       determine initial values for
!      |                                     accuracy check
!      |
!      |
!      Begin parallel region over threads                                 40.22
!      |
!      |-------+ SWOMPU
!      |       |
!      |       |--------+ SWAPAR             determ. of waveparameters
!      |       |
!      |       |--------+ SPROXY             comp. of propagation
!      |       |                             velocities of energy:
!      |       |                             CAX, CAY
!      |       |
!      |       |--------+ SPROSD             comp. of propagation
!      |       |                             velocities of energy:
!      |       |                             CAS, CAD
!      |       |
!      |       |--------+ WINDP1             compute absolute wind, FPM
!      |       |                             mean wind direction, min. and
!      |       |                             max. counters for the wind,
!      |       |                             wind friction velocity
!      |       |
!      |       |--------+ CNTAIL             Compute contributions to
!      |       |                             spectrum due to high frequency
!      |       |                             tail
!      |       |
!      |       |--------+ SPREDT             predict energy density in
!      |       |                             gridpoints for first
!      |       |                             iteration
!      |       |
!      |       |--------+ SINTGRL            comp. Ub, Etot, Hmax, Qb,    40.02
!      |       |        |                    SME, SMA, SMESPC , SMASPC
!      |       |        |
!      |       |        +-------+ FRABRE     comp. of fraction of         30.77
!      |       |                             breaking waves
!      |       |
!      |       |--------+ SOURCE             comp. of source terms
!      |       |        |
!      |       |        +-------+ SBOT       bottom friction
!      |       |        |
!      |       |        +-------+ SVEG       dissipation due to vegetation
!      |       |        |
!      |       |        +-------+ SWCAP      white capping
!      |       |        |
!      |       |        +-------+ SSURF      wave breaking
!      |       |        |
!      |       |        +-------+ SWLTA      nonlinear triad interactions based on LTA
!      |       |        |
!      |       |        +-------+ SWSNL?     nonlinear quadruplet interactions
!      |       |        |
!      |       |        +-------+ SWIND1     first generation wind model
!      |       |                |
!      |       |                + -- WINDP2  compute total wind sea energy o
!      |       |                |    SWIND2  second generation wind model
!      |       |                |
!      |       |                + SWIND3     third generation wind model
!      |       |
!      |       |--------+ ACTION             comp. of ACTION balance eq.
!      |       |        |                    (  @ CAn/@n )
!      |       |        |
!      |       |        +-------+ STIME      comp. of (@AC2/@t)
!      |       |        |
!      |       |        +-------+ STRSX      @[CAX AC2]/@X
!      |       |        |
!      |       |        +-------+ STRSY      @[CAY AC2]/@Y
!      |       |        |
!      |       |        +-------+ STRSS      @[CAS AC2]/@S
!      |       |        |
!      |       |        +-------+ STRSD      @[CAD AC2]/@D
!      |       |
!      |       |--------+ SOLMAT             solve the matrix which is
!      |       |                             filled in SOURCE and ACTION
!      |       |
!      |       |--------+ FILIMP             filter the frequency spectrum
!      |       |   |                         in presence of a current using
!      |       |   |                         a diffusion model (important for
!      |       |   |                         wave blocking) -->IMPLICIT SCHEME
!      |       |   |
!      |       |   |----+ DIFSOL             The matrix filled in FILIMP is
!      |       |                             solved for each direction separately
!      |       |
!      |       |--------+ WINDP3             Limit the energy spectrum
!      |                                     for first and second
!      |                                     generation wind model
!      |
!      |
!      End parallel region over threads                                   40.22
!      |
!      |-------+ PLTSRC                      write sourceterm after an
!      |                                     iteration to a file SOURCE
!      |
!      |-------+ SACCUR / SWSTPC             check accuracy of the comp.
!      |
!     END SWCOMP
!
!
! 12. Structure
!
!     The numerical procedure in SWAN is based on the four-direction
!     point Gauss-Seidel technique with the following
!
!       {**************************************************************}
!       {           definition of the sweep directions                 }
!       {                                                              }
!       {                   \         N         /                      }
!       {     swp_NW = 4     _\|      |      |/_   swp_EN = 3          }
!       {                             |                                }
!       {                             |                                }
!       {             W --------------------------- E  (0 degrees,id=1)}
!       {                             |                                }
!       {                   __.       |      .__                       }
!       {     swp_WS = 1     /|       |      |\    swp_SE = 2          }
!       {                  /          S         \                      }
!       {                                                              }
!       {**************************************************************}
!       {                                                              }
!       { swp_NW:         *  ksy=+1   swp_EN:        *  ksy=+1         }
!       {                 |                          |                 }
!       {                 |  -dy                -dy  |                 }
!       {            dx   |                          |  -dx            }
!       {       *---------o  IX,IY            IX,IY  o--------*        }
!       {     ksx=-1                                    ksx=+1         }
!       {                                                              }
!       {                                                              }
!       { swp_WS:    dx               swp_SE:           -dx            }
!       {       *---------o  IX,IY            IX,IY  o--------*        }
!       {       ksx=-1    |                          |     ksx=+1      }
!       {                 |  dy                  dy  |                 }
!       {                 |                          |                 }
!       {                 *  ksy=-1                  *  ksy=-1         }
!       {**************************************************************}
!
!     ----------------------------------------------------------
!     Call INSAC to give values to HSACC and SACC meant for accuracy check
!     ----------------------------------------------------------
!     For IT = 1 to end of computation time (MTC), do,
!
!       If accuracy <= given tolerance, then do iteration,
!
!         -----------------------------------------------------
!         give argument for sweep : swpdir = 1
!         KSX = -1         DDX = +DX.
!         KSY = -1         DDY = +DY.
!         give number of direction a start and an end value:
!         For IY=2 to MYC and IX=2 to MXC, do,
!            Call SWOMPU to compute the wave field
!         -----------------------------------------------------
!         give argument for sweep : swpdir = 2
!         KSX = +1         DDX = -DX.
!         KSY = -1         DDY = +DY.
!         give number of direction a start and an end value:
!         For IX=MXC-1 to 1 and IY=2 to MYC, do,
!            Call SWOMPU to compute the wave field
!         -----------------------------------------------------
!         give argument for sweep : swpdir = 3
!         KSX = +1.         DDX = -DX.
!         KSY = +1.         DDY = -DY.
!         give number of direction a start and an end value:
!         For IY=MYC-1 to 1 and IX=MXC-1 to 1, do,
!            Call SWOMPU to compute the wave field
!         -----------------------------------------------------
!         give argument for sweep : swpdir = 4
!         KSX = -1         DDX = +DX.
!         KSY = +1         DDY = -DY.
!         give number of direction a start and an end value:
!         For IX=2 to MXC and IY=MYC-1 to 1, do,
!            Call SWOMPU to compute the wave field
!         ----------------------------------------------------
!     CALL PLTSRC to write the source term to a file
!     ----------------------------------------------------
!     CALL SACCUR / SWSTPC to check the accuracy of the computation
!     --------------------------------------------------------
!     End of SWCOMP
!     --------------------------------------------------------
!
! 13. Source text
!
!     ************************************************************************
!     *                                                                      *
!     *                  MAIN SUBROUTINE OF COMPUTATIONAL PART               *
!     *                                                                      *
!     *                               -- SWCOMP --                           *
!     *                                                                      *
!     *                Definition of variables in main program               *
!     *                                                                      *
!     ************************************************************************
!
      INTEGER :: ITER  ,IX    ,IY    ,IS    ,IT                           30.72
      INTEGER :: IP, IDC, ISC                                             40.31
      INTEGER :: KSX   ,KSY   ,SWPDIR
      INTEGER :: INOCNV                                                   30.72
      INTEGER :: INOCNT                                                   40.31
!
      REAL ::  DDX   ,DDY   ,ACCUR ,XIS   ,SNLC1 ,DAL1  ,DAL2  ,DAL3      30.74
!
      LOGICAL :: PRECOR
!
      INTEGER :: MNISL, MXNFL, MXNFR, NPFL, NPFR, NWETP, IDEBUG           40.23
      REAL    :: FRAC                                                     40.23
      PARAMETER (IDEBUG=0)                                                40.23
!
      INTEGER   ISTAT, IF1, IL1                                           40.41
      CHARACTER*20 NUMSTR, CHARS(1)                                       40.41
      CHARACTER*80 MSGSTR                                                 40.41
!
      INTEGER IARR(10)                                                    40.30
      REAL     ARR(10)                                                    40.30

      INTEGER JDUM, JS, JE, JWFRS, JWFRE, INCJ, JJ, JNODE,                40.30
     &        JSD, JED, IE, INCI, II, III, LSTCP                          40.31

      INTEGER IX1, IX2, IY1, IY2                                          40.31
!
!     Add variables for the XNL interface (quadruplet interaction)        40.41
      INTEGER :: IXGRID, IXQUAD, IQERR                                    40.41
!
      INTEGER :: XYTST(2*NPTST)                                           30.21
      INTEGER :: KGRPNT(MXC,MYC)                                          30.21
      INTEGER :: CROSS(2,MCGRD)
!
      REAL     AC2(MDC,MSC,MCGRD)     ,                                   30.21
     &         AC1(MDC,MSC,MCGRD)     ,                                   30.21
     &         COMPDA(MCGRD,MCMVAR)                                       40.31

      REAL WWAWG(8), WWSWG(8)                                             40.22

      INTEGER, DIMENSION(:), ALLOCATABLE :: IDCMIN, IDCMAX,               40.22
     &                                      ISCMIN, ISCMAX                40.22
      INTEGER WWINT(24)                                                   40.41 40.22

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: CAX,CAY,CAX1,CAY1,           40.41 40.22
     &                                       CAS,CAD                      40.22

      REAL, DIMENSION(:,:), ALLOCATABLE :: CGO,KWAVE                      40.22

      REAL, DIMENSION(:,:), ALLOCATABLE :: ALIMW                          40.22

      REAL, DIMENSION(:,:), ALLOCATABLE :: UE,SA1,SA2,SFNL                40.22

      REAL, DIMENSION(:,:), ALLOCATABLE :: DA1C,DA1P,DA1M,                40.22
     &                                     DA2C,DA2P,DA2M,DSNL            40.22

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: MEMNL4                       40.22

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: OBREDF                       40.22
      REAL, DIMENSION(:,:), ALLOCATABLE :: REFLSO                         40.41

      REAL, DIMENSION(:), ALLOCATABLE :: HSAC1,HSAC2,SACC1,SACC2          40.22
      REAL, DIMENSION(:), ALLOCATABLE :: HSAC0,HSDIFC                     40.41
      REAL, DIMENSION(:), ALLOCATABLE :: SACC0,TMDIFC                     40.93

      REAL, DIMENSION(:,:), ALLOCATABLE :: SETPDA                         40.41 40.22

      LOGICAL, DIMENSION(:), ALLOCATABLE :: GROWW                         40.22

      LOGICAL, DIMENSION(:), ALLOCATABLE :: ANYWND                        40.22

      REAL, ALLOCATABLE :: SWTSDA(:,:,:,:)                                40.31

!     SWMATR and LSWMAT replace the single array SWMATR                   40.22
!     that is equivalenced to the logical array LSWMATR                   40.22
!     in the subroutine SWOMPU.                                           40.22
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SWMATR                       40.22

      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LSWMAT                    40.22

!$    LOGICAL, ALLOCATABLE :: LLOCK(:,:)                                  40.31 40.22

      INTEGER, ALLOCATABLE :: ISLMIN(:), NFLIM(:), NRSCAL(:)              40.23

      REAL, ALLOCATABLE :: AC2LOC(:)                                      40.30

!     Add variables for OMP thread parameters.                            40.31
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM        40.31
      INTEGER I1GRD,I2GRD,I1MYC,I2MYC                                     40.31
!
!-----------------------------------------------------------------------
!                      End of variable definition
!-----------------------------------------------------------------------
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOMP')
!
      IF (IT .EQ. 1 .AND. ITEST.GE.1) CALL SWPRSET (SPCSIG)               40.80
!
!     *** print test points ***
!
      IF (NPTST.GT.0) THEN
        DO 121 II = 1, NPTST
          WRITE(PRINTF,1001) II, XYTST(2*II-1)+MXF-2, XYTST(2*II)+MYF-2   40.30
 1001       FORMAT(' Test points :',3I5)
 121    CONTINUE
      ENDIF
!
!     *** prepare ranges of spectral space, constants and      ***        40.41
!     *** weight factors for nonlinear 4 wave interactions     ***        40.41
!                                                                         40.41
!TIMG      CALL SWTSTA(135)                                                    40.41
      IF ( IQUAD.GE.1 )                                                   40.41
     &   CALL FAC4WW (XIS   ,SNLC1 ,                                      40.41
     &                DAL1  ,DAL2  ,DAL3  ,SPCSIG,                        40.41
     &                WWINT ,WWAWG ,WWSWG )                               40.41
!TIMG      CALL SWTSTO(135)                                                    40.41
!                                                                         40.22
! *** Indexing and bounds for SWMAT arrays                                40.22
!                                                                         40.22
      JMATD = 1                                                           40.22
      JMATR = 2                                                           40.22
      JMATL = 3                                                           40.22
      JMATU = 4                                                           40.22
      JMAT5 = 5                                                           40.22
      JMAT6 = 6                                                           40.22
      JDIS0 = 7                                                           40.22
      JDIS1 = JDIS0+MDISP                                                 40.67 40.22
      JGEN0 = JDIS1+MDISP                                                 40.85
      JGEN1 = JGEN0+MGENR                                                 40.85
      JRED0 = JGEN1+MGENR                                                 40.85
      JRED1 = JRED0+MREDS                                                 40.85
      JTRA0 = JRED1+MREDS                                                 40.85
      JTRA1 = JTRA0+MTRNP                                                 40.85
      JAOLD = JTRA1+MTRNP                                                 40.85 40.67 40.22
      JLEK1 = JAOLD+1                                                     40.67 40.22
      MSWMATR = JLEK1                                                     40.85 40.67 40.61 40.22
      JABIN = 1                                                           40.22
      JABLK = 2                                                           40.22
      MLSWMAT = 2                                                         40.22

! *** Stencil size                                                        40.22

      IF (PROPSC.EQ.3) THEN                                               33.09
        ICMAX  = 10                                                       33.08
        LSTCP  = 3                                                        40.31
      ELSE IF (PROPSC.EQ.2) THEN                                          33.10
        ICMAX  = 5                                                        33.10
        LSTCP  = 2                                                        40.31
      ELSE
        ICMAX  = 3                                                        33.10
        LSTCP  = 1                                                        40.31
      ENDIF                                                               33.08

!TIMG      CALL SWTSTA(101)                                                    40.23

!----------------------------------------------------------------------   40.22
!     Begin allocate shared arrays.                                       40.22
!----------------------------------------------------------------------   40.22
!
      ALLOCATE(HSAC1(MCGRD))                                              40.22
      ALLOCATE(HSAC2(MCGRD))                                              40.22
      ALLOCATE(SACC1(MCGRD))                                              40.22
      ALLOCATE(SACC2(MCGRD))                                              40.22
      ALLOCATE(HSAC0(MCGRD))                                              40.41
      ALLOCATE(HSDIFC(MCGRD))                                             40.41
      ALLOCATE(SACC0(MCGRD))                                              40.93
      ALLOCATE(TMDIFC(MCGRD))                                             40.93
!
      ALLOCATE(ISLMIN(MCGRD))                                             40.23
      ALLOCATE(NFLIM(MCGRD))                                              40.23
      ALLOCATE(NRSCAL(MCGRD))                                             40.23
!
      IF ( IQUAD .GE. 1) THEN
!       *** quadruplets ***
        IF ( IQUAD .GE. 3 ) THEN                                          40.10
!         *** prior to every iteration full directional domain ***
          ALLOCATE(MEMNL4(MDC,MSC,MCGRD),STAT=ISTAT)                      40.41 40.22
          IF ( ISTAT.NE.0 ) THEN                                          40.41
             CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')                         40.41
             CALL TXPBLA(CHARS(1),IF1,IL1)                                40.41
             MSGSTR =                                                     40.41
     &         'Allocation problem: array MEMNL4 and return code is '//   40.41
     &         CHARS(1)(IF1:IL1)                                          40.41
             CALL MSGERR ( 4, MSGSTR )                                    40.41
             RETURN                                                       40.41
          END IF                                                          40.41
        ELSE                                                              40.22
!       *** iquad < 3 ***                                                 40.22
          ALLOCATE(MEMNL4(0,0,0))                                         40.22
        END IF
      ELSE
!       *** no quadruplets ***
        if( ALLOCATED(MEMNL4) ) deallocate( MEMNL4)
        ALLOCATE(MEMNL4(0,0,0))                                           40.22
      ENDIF
!
!     *** Lock array for thread management                                40.22
!$    ALLOCATE(LLOCK(MXC,MYC))                                            40.31 40.22
      ALLOCATE(AC2LOC(MCGRD))                                             40.30

      ALLOCATE(SWTSDA(MDC,MSC,NPTSTA,MTSVAR))                             40.31
      SWTSDA = 0.                                                         40.41
!
!     *** In case of SETUP expand array for setup data ***                32.02
!
      IF (LSETUP.GT.0) THEN                                               31.04
        MSTPDA = 23                                                       31.04
        ALLOCATE(SETPDA(MCGRD,MSTPDA))                                    40.41 40.22
      ELSE                                                                40.22
        ALLOCATE(SETPDA(0,0))                                             40.41 40.22
      END IF

!----------------------------------------------------------------------   40.22
!     End allocate shared arrays.                                         40.22
!----------------------------------------------------------------------   40.22

!----------------------------------------------------------------------   40.31
!     Begin initialization of shared arrays.                              40.31
!----------------------------------------------------------------------   40.31

      HSAC1 = 0.                                                          40.31
      HSAC2 = 0.                                                          40.31
      SACC1 = 0.                                                          40.31
      SACC2 = 0.                                                          40.31
      HSAC0 = 0.                                                          40.41
      HSDIFC= 0.                                                          40.41
      SACC0 = 0.                                                          40.93
      TMDIFC= 0.                                                          40.93
      IF ( IQUAD.GE.3 ) MEMNL4 = 0.                                       40.41 40.31
      IF ( LSETUP.GT.0 ) SETPDA = 0.                                      40.31

!----------------------------------------------------------------------   40.31
!     End initialization shared arrays.                                   40.31
!----------------------------------------------------------------------   40.31
!TIMG
!TIMG      CALL SWTSTO(101)                                                    40.23

!----------------------------------------------------------------------   40.22
!     Begin parallel region.                                              40.22
!----------------------------------------------------------------------   40.22

!$OMP PARALLEL DEFAULT(SHARED)                                            40.22
!$OMP+PRIVATE(ITER, SWPDIR, IX, IY, II, IJ, IK)                           40.41 40.22
!$OMP+PRIVATE(CAX, CAY, CAX1, CAY1, CAS, CAD, CGO, KWAVE)                 40.22
!$OMP+PRIVATE(SWMATR, LSWMAT, ALIMW, GROWW, IDCMIN, IDCMAX)               40.22
!$OMP+PRIVATE(ISCMIN, ISCMAX, UE, SA1, SA2, SFNL)                         40.41 40.22
!$OMP+PRIVATE(DA1C, DA1P, DA1M, DA2C, DA2P, DA2M, DSNL)                   40.22
!$OMP+PRIVATE(ANYWND, OBREDF)                                             40.41 40.22
!$OMP+PRIVATE(REFLSO, INOCNT)                                             40.41 40.31
!$OMP+PRIVATE(IP,IDC,ISC)                                                 40.31
!$OMP+PRIVATE(I1GRD,I2GRD,I1MYC,I2MYC)                                    40.31
!$OMP+PRIVATE(JDUM,JJ,III)                                                40.31
!$OMP+PRIVATE(IS,IE,INCI,JS,JE,INCJ,JSD,JED,JNODE,JWFRS,JWFRE)            40.31
!$OMP+COPYIN(ICMAX,CSETUP)
!$OMP+COPYIN(COSLAT,PROPSL)
!$OMP+COPYIN(IPTST,TESTFL)
!
!$OMP MASTER                                                              40.22
!  Print number of threads set by environment                             40.22
!$    IF ( IT.EQ.1 )                                                      41.10
!$   +WRITE(SCREEN,'(a,i2/)')                                             41.10 40.22
!$   +      ' Number of threads during execution of parallel region = ',  41.10 40.22
!$   +      OMP_GET_NUM_THREADS()                                         41.10 40.22
!$OMP END MASTER                                                          40.22
!TIMG
!TIMG      CALL SWTSTA(101)                                                    40.23

!----------------------------------------------------------------------   40.22
!     Begin allocate private arrays.                                      40.22
!----------------------------------------------------------------------   40.22

      ALLOCATE(CAX(MDC,MSC,MICMAX))                                       40.22
      ALLOCATE(CAY(MDC,MSC,MICMAX))                                       40.22
      ALLOCATE(CAX1(MDC,MSC,MICMAX))                                      40.22
      ALLOCATE(CAY1(MDC,MSC,MICMAX))                                      40.22
      ALLOCATE(CAS(MDC,MSC,MICMAX))                                       40.22
      ALLOCATE(CAD(MDC,MSC,MICMAX))                                       40.22
      ALLOCATE(CGO(MSC,MICMAX))                                           40.22
      ALLOCATE(KWAVE(MSC,MICMAX))                                         40.22
!     Since SWMATR has been broken up into a real array(SWMATR) and a     40.22
!     logical array(LSWMAT), the size of each array has been adjusted     40.22
!     to MSWMATR(x-2) and MLSWMAT(2) instead of the original equivalenced 40.22
!     array with a size of MSWMAT(x).                                     40.22
      ALLOCATE(SWMATR(MDC,MSC,MSWMATR))                                   40.22
      ALLOCATE(LSWMAT(MDC,MSC,MLSWMAT))                                   40.22
      ALLOCATE(ALIMW(MDC,MSC))                                            40.22
      ALLOCATE(GROWW(MDC*MSC))                                            40.22
      ALLOCATE(IDCMIN(MSC))                                               40.22
      ALLOCATE(IDCMAX(MSC))                                               40.22
      ALLOCATE(ISCMIN(MDC))                                               40.22
      ALLOCATE(ISCMAX(MDC))                                               40.22
!     *** quadruplets ***
      IF ( IQUAD .GE. 1) THEN
        ALLOCATE(UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                         40.22
        ALLOCATE(SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                        40.22
        ALLOCATE(SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                        40.22
        ALLOCATE(SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                       40.22
        IF ( IQUAD .EQ. 1 ) THEN
!         *** semi-implicit calculation ***
          ALLOCATE(DA1C(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DA1P(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DA1M(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DA2C(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DA2P(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DA2M(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
          ALLOCATE(DSNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))                     40.22
        ELSE                                                              40.22
!       *** iquad > 1 ***                                                 40.22
          ALLOCATE(DA1C(0,0))                                             40.22
          ALLOCATE(DA1P(0,0))                                             40.22
          ALLOCATE(DA1M(0,0))                                             40.22
          ALLOCATE(DA2C(0,0))                                             40.22
          ALLOCATE(DA2P(0,0))                                             40.22
          ALLOCATE(DA2M(0,0))                                             40.22
          ALLOCATE(DSNL(0,0))                                             40.22
        END IF                                                            40.22
      ELSE
!       *** no quadruplets ***
        ALLOCATE(UE(0,0))                                                 40.22
        ALLOCATE(SA1(0,0))                                                40.22
        ALLOCATE(SA2(0,0))                                                40.22
        ALLOCATE(SFNL(0,0))                                               40.22
        ALLOCATE(DA1C(0,0))                                               40.22
        ALLOCATE(DA1P(0,0))                                               40.22
        ALLOCATE(DA1M(0,0))                                               40.22
        ALLOCATE(DA2C(0,0))                                               40.22
        ALLOCATE(DA2P(0,0))                                               40.22
        ALLOCATE(DA2M(0,0))                                               40.22
        ALLOCATE(DSNL(0,0))                                               40.22
      END IF
!
!     *** for wind, indicating a bin inside the wind region     ***
!
      ALLOCATE(ANYWND(MDC))                                               40.22
!
!     *** for obstacles, to store the transmission coefficients ***
!     *** and contribution to the source terms                  ***
!
      ALLOCATE(OBREDF(MDC,MSC,2))                                         40.22
      ALLOCATE(REFLSO(MDC,MSC))                                           40.41

!----------------------------------------------------------------------   40.22
!     End allocate private arrays.                                        40.22
!----------------------------------------------------------------------   40.22

!----------------------------------------------------------------------   40.31
!     Begin initialization of private arrays.                             40.31
!----------------------------------------------------------------------   40.31

      CAX    = 0.                                                         40.31
      CAY    = 0.                                                         40.31
      CAX1   = 0.                                                         40.31
      CAY1   = 0.                                                         40.31
      CAS    = 0.                                                         40.31
      CAD    = 0.                                                         40.31
      CGO    = 0.                                                         40.31
      KWAVE  = 0.                                                         40.31
      SWMATR = 0.                                                         40.31
      ALIMW  = 0.                                                         40.31
      IF ( IQUAD.GE.1 ) THEN                                              40.31
        UE   = 0.                                                         40.31
        SA1  = 0.                                                         40.31
        SA2  = 0.                                                         40.31
        SFNL = 0.                                                         40.31
        IF ( IQUAD.EQ.1 ) THEN                                            40.31
          DA1C = 0.                                                       40.31
          DA1P = 0.                                                       40.31
          DA1M = 0.                                                       40.31
          DA2C = 0.                                                       40.31
          DA2P = 0.                                                       40.31
          DA2M = 0.                                                       40.31
          DSNL = 0.                                                       40.31
        END IF                                                            40.31
      END IF                                                              40.31

!----------------------------------------------------------------------   40.31
!     End initialization private arrays.                                  40.31
!----------------------------------------------------------------------   40.31
!$OMP BARRIER                                                             40.31
!TIMG
!TIMG      CALL SWTSTO(101)                                                    40.23

! Each thread compute its own spatial grid loop bounds for MCGRD          40.31
      CALL SWMTLB(1,MCGRD,I1GRD,I2GRD)                                    40.31
! Each thread compute its own spatial grid loop bounds for MYC            40.31
      CALL SWMTLB(1,MYC,I1MYC,I2MYC)                                      40.31
!
!     *** initialise values for determining the accuracy that ***
!     *** has been reached                                    ***
!     *** This is done in parallel within OpenMP environment  ***         40.31
!
!TIMG      CALL SWTSTA(102)                                                    40.23
      CALL INSAC (AC2               ,SPCSIG          ,COMPDA(1,JDP2)  ,   40.22
     &            HSAC2             ,SACC2           ,KGRPNT          ,   40.31 40.30 40.22
     &            I1MYC             ,I2MYC                            )   40.31
!TIMG      CALL SWTSTO(102)                                                    40.23
!
!     *** To obtain a first estimate of energy density in a    ***
!     *** gridpoint considered we run the SWAN model (in case  ***
!     *** of active wind) in a second generation mode first.   ***
!     *** After 1st iteration, the options, as defined by the  ***
!     *** user, are re-activated.                              ***
!     *** This first guess is not used in nonstationary        ***
!     *** computations (NSTATC>0), or if a restart file was    ***
!     *** used (ICOND=4)                                       ***
!
!$OMP MASTER                                                              40.31
      IF ( IWIND.GE.3 .AND. NSTATC.EQ.0 .AND. ICOND.NE.4 ) THEN           40.41
!     --- first guess will be used
         PRECOR = .TRUE.
      ELSE
         PRECOR = .FALSE.
      END IF
!$OMP END MASTER                                                          40.31
!
!     *** call initialization procedure of XNL to create *.BQF ***        40.41
!     *** interaction files                                    ***        40.41
!                                                                         40.41
!TIMG      CALL SWTSTA(135)                                                    40.41
!$OMP MASTER
      IF (IQUAD.EQ.51.OR.IQUAD.EQ.52.OR.IQUAD.EQ.53) THEN                 40.41
         CALL init_constants                                              40.41
         IXQUAD = IQUAD - 50                                              40.41
         IF (IAMMASTER) THEN                                              40.95 40.41
            WRITE(SCREEN,*) 'GurboQuad initialization'                    40.41
            WRITE(SCREEN,*) 'gravity               :', GRAV               40.41
            WRITE(SCREEN,*) 'pftail                :', PWTAIL(1)          40.41
            WRITE(SCREEN,*) 'number of sigma values:', MSC                40.41
            WRITE(SCREEN,*) 'number of directions  :', MDC                40.41
            WRITE(SCREEN,*) 'IQ_QUAD               :', IXQUAD             40.41
         END IF                                                           40.41
!                                                                         40.41
         IXGRID = 3                                                       40.41
         CALL xnl_init(SPCSIG  , SPCDIR(:,1)*180./PI , MSC , MDC ,        40.41
     &                 -PWTAIL(1), GRAV   , COMPDA(2,JDP2) ,              40.41
     &                 MCGRD-1   , IXQUAD , IXGRID   ,INODE,IQERR )       40.41
      END IF                                                              40.41
!$OMP END MASTER
!TIMG      CALL SWTSTO(135)                                                    40.41
!

!TIMG      CALL SWTSTA(103)                                                    40.23
      DO 450 ITER = 1, ITERMX                                             30.00

!       initialise local (thread private) counter for SIP solver          40.31
        INOCNT = 0                                                        40.31
!
!       initialise propagation, generation, dissipation, redistribution,  40.85
!       leak and radiation stress for each iteration
!       this is done in parallel within OpenMP environment                40.31
!
        DO IP = I1GRD,I2GRD                                               40.31
          COMPDA(IP,JDISS) = 0.
          COMPDA(IP,JLEAK) = 0.
          COMPDA(IP,JDSXB) = 0.                                           40.61
          COMPDA(IP,JDSXS) = 0.                                           40.61
          COMPDA(IP,JDSXW) = 0.                                           40.61
          COMPDA(IP,JDSXV) = 0.                                           40.61
          COMPDA(IP,JGENR) = 0.                                           40.85
          COMPDA(IP,JGSXW) = 0.                                           40.85
          COMPDA(IP,JREDS) = 0.                                           40.85
          COMPDA(IP,JRSXQ) = 0.                                           40.85
          COMPDA(IP,JRSXT) = 0.                                           40.85
          COMPDA(IP,JTRAN) = 0.                                           40.85
          COMPDA(IP,JTSXG) = 0.                                           40.85
          COMPDA(IP,JTSXT) = 0.                                           40.85
          COMPDA(IP,JTSXS) = 0.                                           40.85
          COMPDA(IP,JRADS) = 0.                                           40.85
          COMPDA(IP,JQB  ) = 0.                                           40.67
        ENDDO
!
!       initialise Ursell number to 0 for each iteration                  40.03
!       this is done in parallel within OpenMP environment                40.31
!
        IF (ITRIAD.GT.0
     &                 ) THEN
           DO IP = I1GRD,I2GRD                                            40.31
            COMPDA(IP,JURSEL) = 0.
          ENDDO
        ENDIF
!                                                                         40.31
!       *** IQUAD = 3: the nonlinear wave interactions are     ***        40.31
!       *** calculated just once for an iteration. First,      ***        40.31
!       *** set the auxiliary array equal zero before a        ***        40.31
!       *** new iteration                                      ***        40.31
!       *** This is done in parallel within OpenMP environment ***        40.31
!                                                                         40.31
        IF ( IQUAD .GE. 3 ) THEN                                          40.31
          DO IP = I1GRD,I2GRD                                             40.31
            DO ISC = 1,MSC                                                40.31
              DO IDC = 1,MDC                                              40.31
                MEMNL4(IDC,ISC,IP)=0.                                     40.31
              END DO                                                      40.31
            END DO                                                        40.31
          END DO                                                          40.31
        END IF                                                            40.31
!                                                                         40.31
!----------------------------------------------------------------------   40.31
!     Begin master thread region.                                         40.31
!----------------------------------------------------------------------   40.31
!$OMP MASTER                                                              40.31
!                                                                         40.31
!       *** If a current is present and a penta-diagonal solver    ***    40.31
!       *** is employed, it is possible that the solver does       ***    40.31
!       *** not converged. For this, the counter INOCNV represents ***    40.31
!       *** the number of geographical points in which the solver  ***    40.31
!       *** did not converged                                      ***    40.31
!                                                                         40.31
        INOCNV = 0                                                        40.31
!
        IF ( ITEST.GE.30 .OR. IDEBUG.EQ.1 ) THEN
           ISLMIN(:) = 9999                                               40.23
           NFLIM (:) = 0                                                  40.23
           NRSCAL(:) = 0                                                  40.23
           IARR = 0                                                       40.30
           ARR  = 0.                                                      40.30
        END IF
!
        IF ( PRECOR ) THEN                                                40.41
!           *** third generation wave input ***
           IF ( ITER .EQ. 1 )THEN
!             *** save settings of 3rd generation model             ***   40.13
!             *** bottom friction, surf breaking and triads may     ***
!             *** still active                                      ***
              KWIND  = IWIND
              KWCAP  = IWCAP
              KQUAD  = IQUAD
!             ***  save maximum change per bin and under-relaxation ***   40.41 40.13
              GRWOLD = PNUMS(20)
              ALFAT  = PNUMS(30)                                          40.31

!             first guess settings                                        40.13

!             if 1st generation is to be used as first guess, replace     40.13
!             the next statement by IWIND = 1                             40.13
              IWIND  = 2
              IWCAP  = 0
              IQUAD  = 0
              PNUMS(20) = 1.E22
!             ***  under-relaxation parameter is PNUMS(30) and            40.31
!                  temporarily set to zero                                40.31
              PNUMS(30) = 0.                                              40.31
           ELSE IF ( ITER .EQ. 2 ) THEN
              IWIND  = KWIND
              IWCAP  = KWCAP
              IQUAD  = KQUAD
              PNUMS(20) = GRWOLD
              PNUMS(30) = ALFAT                                           40.31
           ENDIF
!
        ENDIF

        IF ( PRECOR .AND. ITER .LE. 2 ) THEN                              40.41
           WRITE(PRINTF,*)' -----------------------------------------',
     &                    '----------------------'
           IF ( ITER .EQ. 1 ) THEN
              WRITE(PRINTF,*) ' First guess by 2nd generation model',
     &                        ' flags for first iteration:'
           ELSE IF ( ITER .EQ. 2 ) THEN
              WRITE(PRINTF,*) ' Options given by user are activated',
     &                        ' for proceeding calculation:'
           ENDIF
           WRITE(PRINTF,2001) ITER, PNUMS(20), PNUMS(30)
 2001      FORMAT('  ITER    ',I4,' GRWMX    ',E12.4,' ALFA     ',
     &            E12.4)
           WRITE(PRINTF,2002) IWIND, IWCAP, IQUAD
 2002      FORMAT('  IWIND   ',I4,' IWCAP   ',I4,' IQUAD   ',I4)
           WRITE(PRINTF,2003) ITRIAD, IBOT , ISURF
 2003      FORMAT('  ITRIAD  ',I4,' IBOT    ',I4,' ISURF   ',I4)
           WRITE(PRINTF,2005) IVEG                                        40.55
 2005      FORMAT('  IVEG    ',I4)                                        40.55
           WRITE(PRINTF,*)' -----------------------------------------',
     &                    '----------------------'
        ENDIF

!       --- calculate diffraction parameter and its derivatives           40.21
!TIMG        CALL SWTSTA(137)                                                  40.41
        IF ( IDIFFR.GT.0 )                                                40.21
     &     CALL DIFPAR( AC2   , SPCSIG, KGRPNT, COMPDA(1,JDP2),           40.21
     &                  CROSS , XCGRID, YCGRID, XYTST  )                  40.21
!TIMG        CALL SWTSTO(137)                                                  40.41

!----------------------------------------------------------------------   40.22
!     End master thread region.                                           40.22
!----------------------------------------------------------------------   40.22
!$OMP END MASTER                                                          40.22

!----------------------------------------------------------------------   40.22
!     Synchronize threads before loop over sweep directions.              40.22
!----------------------------------------------------------------------   40.22
!$OMP BARRIER                                                             40.22
!
!               *** START ITERATION PROCESS WITH 4 SWEEPS ***
!
!       *** loop over sweep directions ***
!
        DO 410 SWPDIR = 1, 4

!           Initialize LLOCK in parallel                                  40.31
!           Make .FALSE. at grid points where depth is negative           40.31
!$          DO IY = I1MYC,I2MYC                                           40.31
!$            DO IX = 1, MXC                                              40.31
!$              IF (COMPDA(KGRPNT(IX,IY),JDP2).GE.DEPMIN) THEN            40.31
!$                 LLOCK(IX,IY) = .TRUE.                                  40.31
!$              ELSE                                                      40.31
!$                 LLOCK(IX,IY) = .FALSE.                                 40.31
!$              ENDIF                                                     40.31
!$            ENDDO                                                       40.31
!$          ENDDO                                                         40.31

!----------------------------------------------------------------------   40.31
!     Synchronize threads before setting LLOCK for boundary.              40.31
!----------------------------------------------------------------------   40.31
!$OMP BARRIER                                                             40.31

!----------------------------------------------------------------------   40.31
!     Begin master thread region.                                         40.31
!----------------------------------------------------------------------   40.31
!$OMP MASTER                                                              40.31

!           make LLOCK False for points on boundary                       40.22
            IF (SWPDIR.EQ.1) THEN
              KSX = -1
              KSY = -1
              DDX = +DX
              DDY = +DY
              IF (KREPTX.EQ.0) THEN
                IX1 = 2
                IF (.NOT.LMXF) IX1 = IX1-1+IHALOX                         40.41 40.31
!$              LLOCK(IX1-1,:) = .FALSE.                                  40.31 40.22
              ELSE
                IX1 = 1
              ENDIF
              IX2 = MXC
              IY1 = 2
              IY2 = MYC
              IF (.NOT.LMXL) IX2 = IX2-IHALOX                             40.41 40.31
              IF (.NOT.LMYF) IY1 = IY1-1+IHALOY                           40.41 40.31
              IF (.NOT.LMYL) IY2 = IY2-IHALOY                             40.41 40.31
!$            LLOCK(:,IY1-1) = .FALSE.                                    40.31 40.22
            ELSE IF (SWPDIR.EQ.2) THEN
              KSX = +1
              KSY = -1
              DDX = -DX
              DDY = +DY
              IF (KREPTX.EQ.0) THEN
                IX1 = MXC-1
                IF (.NOT.LMXL) IX1 = IX1+1-IHALOX                         40.41 40.31
!$              LLOCK(IX1+1,:) = .FALSE.                                  40.31 40.22
              ELSE
                IX1 = MXC
              ENDIF
              IX2 = 1
              IY1 = 2
              IY2 = MYC
              IF (.NOT.LMXF) IX2 = IX2+IHALOX                             40.41 40.31
              IF (.NOT.LMYF) IY1 = IY1-1+IHALOY                           40.41 40.31
              IF (.NOT.LMYL) IY2 = IY2-IHALOY                             40.41 40.31
!$            LLOCK(:,IY1-1) = .FALSE.                                    40.31 40.22
            ELSE IF (SWPDIR.EQ.3) THEN
              KSX = +1
              KSY = +1
              DDX = -DX
              DDY = -DY
              IF (KREPTX.EQ.0) THEN
                IX1 = MXC-1
                IF (.NOT.LMXL) IX1 = IX1+1-IHALOX                         40.41 40.31
!$              LLOCK(IX1+1,:) = .FALSE.                                  40.31 40.22
              ELSE
                IX1 = MXC
              ENDIF
              IX2 = 1
              IY1 = MYC-1
              IY2 = 1
              IF (.NOT.LMXF) IX2 = IX2+IHALOX                             40.41 40.31
              IF (.NOT.LMYL) IY1 = IY1+1-IHALOY                           40.41 40.31
              IF (.NOT.LMYF) IY2 = IY2+IHALOY                             40.41 40.31
!$            LLOCK(:,IY1+1) = .FALSE.                                    40.31 40.22
            ELSE IF (SWPDIR.EQ.4) THEN
              KSX = -1
              KSY = +1
              DDX = +DX
              DDY = -DY
              IF (KREPTX.EQ.0) THEN
                IX1 = 2
                IF (.NOT.LMXF) IX1 = IX1-1+IHALOX                         40.41 40.31
!$              LLOCK(IX1-1,:) = .FALSE.                                  40.31 40.22
              ELSE
                IX1 = 1
              ENDIF
              IX2 = MXC
              IY1 = MYC-1
              IY2 = 1
              IF (.NOT.LMXL) IX2 = IX2-IHALOX                             40.41 40.31
              IF (.NOT.LMYL) IY1 = IY1+1-IHALOY                           40.41 40.31
              IF (.NOT.LMYF) IY2 = IY2+IHALOY                             40.41 40.31
!$            LLOCK(:,IY1+1) = .FALSE.                                    40.31 40.22
            ENDIF
!
            IYSTEP = KSY
!
!           *** change values of variables for one-dimensional run ***    32.02
!
            IF ( ONED ) THEN                                              32.02
              IY1    = 1                                                  32.02
              IY2    = 1                                                  32.02
              KSY    = 0                                                  32.02
              IYSTEP = 1
            ENDIF                                                         32.02
!
            IF (SCREEN.NE.PRINTF) THEN
              IF (NSTATC.EQ.1) THEN                                       40.00
!                IF (IAMMASTER) WRITE(SCREEN,313) CHTIME, IT,              40.95 40.30
!     &                                                 ITER, SWPDIR       40.30
 313            FORMAT ('+time ', A18, ', step ',I6, '; iteration '
     &                  ,I4, '; sweep ',I1)                               40.23 40.00
              ELSE
                WRITE(PRINTF,314) ITER, SWPDIR
                IF (IAMMASTER) THEN                                       40.95 40.30
                   IF (SWPDIR.EQ.1) THEN                                  40.30
!                      WRITE(SCREEN,314) ITER, SWPDIR                      40.30
                   ELSE                                                   40.30
!                      WRITE(SCREEN,315) ITER, SWPDIR                      40.30
                   END IF                                                 40.30
                END IF                                                    40.30
 314            FORMAT (' iteration ', I4, '; sweep ', I1)                40.23 30.50
 315            FORMAT ('+iteration ', I4, '; sweep ', I1)                40.23 30.50
              ENDIF
            ENDIF
!
!----------------------------------------------------------------------   40.31
!     End master thread region and synchronize threads                    40.31
!----------------------------------------------------------------------   40.31
!$OMP END MASTER                                                          40.31
!$OMP BARRIER                                                             40.31
!$OMP FLUSH                                                               40.31
!
! ======================================================================  40.30
!                                                                         40.30
!           Set up start and end indices for respectively                 40.30
!           IX- and IY-loops appropriated for block                       40.30
!           wavefront approach within distributed-memory                  40.30
!           environment                                                   40.30
!                                                                         40.30
! ======================================================================  40.30

            IF ( MXCGL.GT.MYCGL .OR. .NOT.PARLL ) THEN                    40.30
               JS   =  IY1                                                40.30
               JE   =  IY2                                                40.30
               INCJ = -IYSTEP                                             40.30
               IS   =  IX1                                                40.30
               IE   =  IX2                                                40.30
               INCI = -KSX                                                40.30
               IF (SWPDIR.EQ.1) THEN                                      40.30
                  JSD   = JE - 1                                          40.30
                  JED   = JS - 1                                          40.30
                  JNODE = INODE                                           40.30
                  JWFRS = 0                                               40.30
                  JWFRE = NPROC-1                                         40.30
               ELSE IF (SWPDIR.EQ.2) THEN                                 40.30
                  JSD   = JE - 1                                          40.30
                  JED   = JS - 1                                          40.30
                  JNODE = NPROC+1-INODE                                   40.30
                  JWFRS = 0                                               40.30
                  JWFRE = NPROC-1                                         40.30
               ELSE IF (SWPDIR.EQ.3) THEN                                 40.30
                  JSD   = JS - 1                                          40.30
                  JED   = JE - 1                                          40.30
                  JNODE = NPROC+1-INODE                                   40.30
                  JWFRS = NPROC-1                                         40.30
                  JWFRE = 0                                               40.30
               ELSE IF (SWPDIR.EQ.4) THEN                                 40.30
                  JSD   = JS - 1                                          40.30
                  JED   = JE - 1                                          40.30
                  JNODE = INODE                                           40.30
                  JWFRS = NPROC-1                                         40.30
                  JWFRE = 0                                               40.30
               END IF                                                     40.30
            ELSE                                                          40.30
               JS   =  IX1                                                40.30
               JE   =  IX2                                                40.30
               INCJ = -KSX                                                40.30
               IS   =  IY1                                                40.30
               IE   =  IY2                                                40.30
               INCI = -IYSTEP                                             40.30
               IF (SWPDIR.EQ.1) THEN                                      40.30
                  JSD   = JE - 1                                          40.30
                  JED   = JS - 1                                          40.30
                  JNODE = INODE                                           40.30
                  JWFRS = 0                                               40.30
                  JWFRE = NPROC-1                                         40.30
               ELSE IF (SWPDIR.EQ.2) THEN                                 40.30
                  JSD   = JS - 1                                          40.30
                  JED   = JE - 1                                          40.30
                  JNODE = INODE                                           40.30
                  JWFRS = NPROC-1                                         40.30
                  JWFRE = 0                                               40.30
               ELSE IF (SWPDIR.EQ.3) THEN                                 40.30
                  JSD   = JS - 1                                          40.30
                  JED   = JE - 1                                          40.30
                  JNODE = NPROC+1-INODE                                   40.30
                  JWFRS = NPROC-1                                         40.30
                  JWFRE = 0                                               40.30
               ELSE IF (SWPDIR.EQ.4) THEN                                 40.30
                  JSD   = JE - 1                                          40.30
                  JED   = JS - 1                                          40.30
                  JNODE = NPROC+1-INODE                                   40.30
                  JWFRS = 0                                               40.30
                  JWFRE = NPROC-1                                         40.30
               END IF                                                     40.30
            END IF                                                        40.30

!----------------------------------------------------------------------   40.22
!     Execute loop over rows of spatial grid in a                         40.22
!     pipelined parallel manner within OpenMP environment                 40.22
!----------------------------------------------------------------------   40.22
!$OMP DO SCHEDULE(STATIC,1)                                               40.22
!$OMP+FIRSTPRIVATE(WWINT)                                                 40.31 40.22
!$OMP+LASTPRIVATE(WWINT)                                                  40.31 40.22

! ======================================================================  40.30
!                                                                         40.30
!           Within distributed-memory environment, current                40.30
!           sweep is carry out in a block wavefront manner                40.30
!                                                                         40.30
! ======================================================================  40.30

            DO 400 JDUM = JS+JWFRS, JE+JWFRE, INCJ                        40.30

             JJ = JDUM-JNODE+1                                            40.30
             IF (JNODE.GE.JDUM-JSD .AND. JNODE.LE.JDUM-JED) THEN          40.30

! ======================================================================  40.30
!                                                                         40.30
!             Receive action density from previous updated                40.30
!             row within distributed-memory environment                   40.30
!                                                                         40.30
! ======================================================================  40.30

!TIMG!MPI              CALL SWTSTA(213)                                            40.30
              IF ( MXCGL.GT.MYCGL ) THEN                                  40.30
                 DO III = LSTCP, 1, -1                                    40.31
                    CALL SWRECVAC(AC2,IS-III*INCI,JJ,SWPDIR,KGRPNT)       40.31
                 END DO                                                   40.31
              ELSE                                                        40.30
                 DO III = LSTCP, 1, -1                                    40.31
                    CALL SWRECVAC(AC2,JJ,IS-III*INCI,SWPDIR,KGRPNT)       40.31
                 END DO                                                   40.31
              END IF                                                      40.30
!TIMG!MPI              CALL SWTSTO(213)                                            40.30
              IF (STPNOW()) RETURN                                        40.30

              DO 390 II = IS, IE, INCI                                    40.30
                IF ( MXCGL.GT.MYCGL .OR. .NOT.PARLL ) THEN                40.30
                  IX = II                                                 40.30
                  IY = JJ                                                 40.30
                ELSE                                                      40.30
                  IX = JJ                                                 40.30
                  IY = II                                                 40.30
                END IF                                                    40.30

!----------------------------------------------------------------------   40.22
!               The next while loop will guarantee execution within       40.22
!               OpenMP environment will not proceed until the data        40.22
!               dependencies for grid point (IX,IY) are satisfied.        40.22
!               Since we parallelize only in the y-direction, we only     40.31
!               need to check data dependencies in the y-direction.       40.31
!               The flush is required to ensure each thread has a         40.22
!               consistent view of LLOCK.                                 40.22
!----------------------------------------------------------------------   40.22
!$              IF ( .NOT.ONED ) THEN                                     40.31
!$                  DO WHILE(LLOCK(IX,IY+IYSTEP))                         40.31
!$OMP FLUSH (LLOCK)                                                       40.31
!$                  END DO                                                40.31
!$OMP FLUSH                                                               40.31
!$              END IF                                                    40.31

!TIMG                CALL SWTSTA(104)                                          40.23
                CALL SWOMPU (SWPDIR,KSX              ,KSY              ,
     &            IX               ,IY               ,DDX              ,
     &            DDY              ,DT               ,SNLC1            ,
     &            DAL1             ,DAL2             ,DAL3             ,
     &            XIS              ,SWTSDA           ,INOCNT           ,  40.31
     &            AC2              ,COMPDA           ,SPCDIR           ,
     &            SPCSIG           ,XYTST            ,ITER             ,  30.72
     &                              CGO              ,                    40.41 40.22
     &            CAX              ,CAY              ,CAS              ,  40.22
     &            CAD              ,SWMATR           ,LSWMAT           ,  40.22
     &            KWAVE            ,                                      40.22
     &            ALIMW            ,GROWW                              ,  40.17 40.22
     &            UE               ,SA1              ,SA2              ,  40.22
     &            DA1C             ,DA1P             ,DA1M             ,  40.22
     &            DA2C             ,DA2P             ,DA2M             ,  40.22
     &            SFNL             ,DSNL             ,MEMNL4           ,  40.22
     &            IDCMIN           ,IDCMAX           ,                    40.41 40.22
     &            WWINT            ,WWAWG            ,WWSWG            ,  40.22
     &            ISCMIN           ,ISCMAX           ,                    40.41 40.22
     &            ANYWND           ,AC1              ,IT               ,  40.22
     &            XCGRID           ,YCGRID           ,                    40.41 30.72
     &            KGRPNT           ,CROSS            ,                    300597
     &            OBREDF           ,REFLSO           ,                    40.41 40.22
     &            ISLMIN           ,NFLIM            ,NRSCAL           ,  40.23
     &            CAX1             ,CAY1                                  40.22
     &                                                                 )
!TIMG                CALL SWTSTO(104)                                          40.23
                IF (STPNOW()) RETURN                                      34.01

!----------------------------------------------------------------------   40.22
!               Once the computation is done for grid point (IX,IY) the   40.22
!               thread signals that the data is available by changing     40.22
!               LLOCK(IX,IY).                                             40.22
!----------------------------------------------------------------------   40.22
!$OMP FLUSH                                                               40.31
!$              LLOCK(IX,IY) = .FALSE.                                    40.31 40.22
!$OMP FLUSH (LLOCK)                                                       40.31

 390          CONTINUE

! ======================================================================  40.30
!                                                                         40.30
!             Send action density to next row                             40.30
!             within distributed-memory environment                       40.30
!                                                                         40.30
! ======================================================================  40.30

!TIMG!MPI              CALL SWTSTA(213)                                            40.30
              IF ( MXCGL.GT.MYCGL ) THEN                                  40.30
                 DO III = LSTCP-1, 0, -1                                  40.31
                    CALL SWSENDAC(AC2,IE-III*INCI,JJ,SWPDIR,KGRPNT)       40.31
                 END DO                                                   40.31
              ELSE                                                        40.30
                 DO III = LSTCP-1, 0, -1                                  40.31
                    CALL SWSENDAC(AC2,JJ,IE-III*INCI,SWPDIR,KGRPNT)       40.31
                 END DO                                                   40.31
              END IF                                                      40.30
!TIMG!MPI              CALL SWTSTO(213)                                            40.30
              IF (STPNOW()) RETURN                                        40.30
             END IF                                                       40.30

 400        CONTINUE
!$OMP ENDDO NOWAIT                                                        40.22

!----------------------------------------------------------------------   40.22
!     Synchronize threads before checking stop condition and              40.22
!     before starting next sweep direction.                               40.22
!----------------------------------------------------------------------   40.22
!$OMP BARRIER                                                             40.22
 410    CONTINUE
!
!       --- exchange action densities at subdomain interfaces             40.31
!           within distributed-memory environment if S&L scheme           40.31
!           is employed because of one downwind point present             40.31
!
        IF (PROPSC.EQ.3) THEN                                             40.31
!TIMG!MPI           CALL SWTSTA(213)                                               40.31
           DO ID = 1, MDC                                                 40.31
              DO IS = 1, MSC                                              40.31
                 AC2LOC(:) = AC2(ID,IS,:)                                 40.31
                 CALL SWEXCHG( AC2LOC, KGRPNT )                           40.31
                 AC2(ID,IS,:) = AC2LOC(:)                                 40.31
              END DO                                                      40.31
           END DO                                                         40.31
!TIMG!MPI           CALL SWTSTO(213)                                               40.31
           IF (STPNOW()) RETURN                                           40.31
        END IF                                                            40.31

!----------------------------------------------------------------------   40.31
!     Each thread sum contributions to the global INOCNV counter          40.31
!     which counts the number of grid points over the four sweeps         40.31
!     in which the SIP solver did not converge.                           40.31
!----------------------------------------------------------------------   40.31
!$OMP ATOMIC                                                              40.31
      INOCNV = INOCNV + INOCNT                                            40.31

!----------------------------------------------------------------------   40.31
!     Synchronize threads before master thread stores source terms        40.31
!     and computes wave induced setup.                                    40.31
!----------------------------------------------------------------------   40.31
!$OMP BARRIER                                                             40.31

!----------------------------------------------------------------------   40.22
!     Begin master thread region.                                         40.22
!----------------------------------------------------------------------   40.22
!$OMP MASTER                                                              40.22
!
        IF ( ITEST.GE.30 .OR. IDEBUG.EQ.1 ) THEN                          40.23
           NWETP = 0                                                      40.23
           DO 109 IP = 2, MCGRD                                           40.23
              IF (COMPDA(IP,JDP2).GT.DEPMIN) NWETP = NWETP + 1            40.23
 109       CONTINUE                                                       40.23
!                                                                         40.23
           MXNFL = MAXVAL(NFLIM)                                          40.23
           NPFL  = COUNT(MASK=NFLIM>0)                                    40.23
!                                                                         40.23
           MNISL = MINVAL(ISLMIN)                                         40.23
!                                                                         40.23
           MXNFR = MAXVAL(NRSCAL)                                         40.23
           NPFR  = COUNT(MASK=NRSCAL>0)                                   40.23
!                                                                         40.23
! ======================================================================  40.30
!                                                                         40.30
!          Gather data meant for global reductions in arrays              40.30
!          IARR and ARR within distributed-memory environment             40.30
!                                                                         40.30
! ======================================================================  40.30
!                                                                         40.30
           IARR(1) = NWETP                                                40.30
           IARR(2) = NPFL                                                 40.30
           IARR(3) = NPFR                                                 40.30
           IARR(4) = INOCNV                                               40.30
!                                                                         40.30
           ARR(1)  = REAL(MXNFL)                                          40.30
           ARR(2)  = REAL(MXNFR)                                          40.30
!                                                                         40.30
!          --- carry out reductions across all nodes                      40.30
!                                                                         40.30
           CALL SWREDUCE(   ARR, 4, SWREAL, SWMAX )                       40.30
           IF (STPNOW()) RETURN                                           40.30
           CALL SWREDUCE(  IARR, 4, SWINT , SWSUM )                       40.30
           IF (STPNOW()) RETURN                                           40.30
           CALL SWREDUCE( MNISL, 1, SWINT , SWMIN )                       40.30
           IF (STPNOW()) RETURN                                           40.30
!                                                                         40.30
           NWETP     = IARR(1)                                            40.30
           NPFL      = IARR(2)                                            40.30
           NPFR      = IARR(3)                                            40.30
           INOCNV    = IARR(4)                                            40.30
!                                                                         40.30
           MXNFL     = NINT(ARR(1))                                       40.30
           MXNFR     = NINT(ARR(2))                                       40.30
!                                                                         40.30
           FRAC = REAL(NPFL)*100./REAL(NWETP)                             40.23
           IF(NPFL.GT.0) WRITE(PRINTF,110) 'limiter',FRAC,MXNFL           40.23
           IF (NSTATC.EQ.0 .AND. NPFL.GT.0 .AND. IAMMASTER)               40.95 40.23
     &        WRITE(SCREEN,110) 'limiter',FRAC,MXNFL                      40.23
!                                                                         40.23
           IF (NPFL.GT.0) WRITE(PRINTF,111) SPCSIG(MNISL)/PI2             40.23
           IF (NSTATC.EQ.0 .AND. NPFL.GT.0 .AND. IAMMASTER)               40.95 40.23
     &        WRITE(SCREEN,111) SPCSIG(MNISL)/PI2                         40.23
!                                                                         40.23
           FRAC = REAL(NPFR)*100./REAL(NWETP)                             40.23
           IF(NPFR.GT.0) WRITE(PRINTF,110) 'rescaling',FRAC,MXNFR         40.23
           IF (NSTATC.EQ.0 .AND. NPFR.GT.0 .AND. IAMMASTER)               40.95 40.23
     &        WRITE(SCREEN,110) 'rescaling',FRAC,MXNFR                    40.23
!                                                                         40.23
 110       FORMAT(1X,'use of ',A9,' in ',F6.2,                            40.23
     &         ' % of wet points with maximum in spectral space = ',      40.23
     &         I4)                                                        40.23
 111       FORMAT(1X,                                                     40.23
     &     'lowest frequency occured above which limiter is applied = ',  40.23
     &     F7.4,' Hz')                                                    40.23
        END IF                                                            40.23
!
!       *** store the source terms for test gridpoints  ***
!       *** in the files IFPAR, IFS1D and IFS2D         ***
!
!TIMG        CALL SWTSTA(105)                                                  40.23
        IF (NPTST.GT.0 .AND. NSTATM.EQ.0                                  40.00
     &                  ) THEN
          IF (IFPAR.GT.0) WRITE (IFPAR, 12) ITER                          40.00
          IF (IFS1D.GT.0) WRITE (IFS1D, 12) ITER                          40.00
          IF (IFS2D.GT.0) WRITE (IFS2D, 12) ITER                          40.00
  12      FORMAT (I4, T41, 'iteration')
          CALL PLTSRC (SWTSDA(1,1,1,JPWNDA)  ,SWTSDA(1,1,1,JPWNDB)  ,
     &                 SWTSDA(1,1,1,JPWCAP)  ,SWTSDA(1,1,1,JPBTFR)  ,
     &                 SWTSDA(1,1,1,JPWBRK)  ,SWTSDA(1,1,1,JP4S)    ,
     &                 SWTSDA(1,1,1,JP4D)    ,SWTSDA(1,1,1,JPTRI)   ,
     &                 SWTSDA(1,1,1,JPVEGT)  ,                            40.55
     &                 AC2                   ,SPCSIG                ,     40.00
     &                 COMPDA(1,JDP2)        ,XYTST                 ,
     &                                        KGRPNT                )     40.00
        END IF
!TIMG        CALL SWTSTO(105)                                                  40.23
!
!       *** compute wave-induced setup ***                                32.02
!
!TIMG        CALL SWTSTA(106)                                                  40.23
        IF (LSETUP.GT.0)                                                  31.03
     &    CALL SETUPP ( KGRPNT, MSTPDA, SETPDA, AC2, COMPDA(1,JDP2),      40.41
     &                  COMPDA(1,JDPSAV), COMPDA(1,JSETUP),               40.41
     &                  XCGRID, YCGRID, SPCSIG, SPCDIR )                  40.41
!TIMG        CALL SWTSTO(106)                                                  40.23
!
!----------------------------------------------------------------------   40.31
!     End master thread region.                                           40.31
!----------------------------------------------------------------------   40.31
!$OMP END MASTER                                                          40.31
!
!       *** check if numerical accuracy has been reached       ***
!       *** this is done in parallel within OpenMP environment ***        40.31
!
!TIMG        CALL SWTSTA(102)                                                  40.23
        IF (PNUMS(21).EQ.0.) THEN                                         40.41
           CALL SACCUR (COMPDA(1,JDP2),KGRPNT          ,                  40.30
     &             XYTST           ,                                      40.41
     &             AC2             ,SPCSIG          ,ACCUR           ,    30.72
     &             HSAC1           ,HSAC2           ,SACC1           ,    30.90 40.22
     &             SACC2           ,COMPDA(1,JDHS)  ,COMPDA(1,JDTM)  ,    40.31 30.90 40.22
     &             I1MYC           ,I2MYC                            )    40.31
        ELSE IF (PNUMS(21).EQ.1.) THEN                                    40.41
           CALL SWSTPC ( HSAC0         ,HSAC1           ,HSAC2 ,          40.41
     &                   SACC0         ,SACC1           ,SACC2 ,          40.93 40.41
     &                   HSDIFC        ,TMDIFC          ,                 40.93 40.41
     &                   COMPDA(1,JDHS),COMPDA(1,JDTM)  ,                 40.41
     &                   COMPDA(1,JDP2),ACCUR           ,                 40.41
     &                   I1MYC         ,I2MYC           )                 40.93
        END IF                                                            40.41
!TIMG        CALL SWTSTO(102)                                                  40.23
!
!----------------------------------------------------------------------   40.31
!     Begin master thread region.                                         40.31
!----------------------------------------------------------------------   40.31
!$OMP MASTER                                                              40.31
!
!       *** info regarding the iteration process and the accuracy ***
!
        IF (PNUMS(21).EQ.1. .AND. ITER.EQ.1) THEN                         40.41
           WRITE(PRINTF,113)                                              40.41
           IF (NSTATC.EQ.0.AND.IAMMASTER) WRITE(SCREEN,113)               40.95 40.41
        ELSE                                                              40.41
           WRITE(PRINTF,112) ACCUR !,PNUMS(4)                               30.72
           IF (NSTATC.EQ.0.AND.IAMMASTER)                                 40.95 40.30
     &        WRITE(SCREEN,112) ACCUR !,PNUMS(4)                            40.30 40.00
        END IF
 112    FORMAT('+',' waves circulation completion: ',F10.3,' %')                                            40.41
 113    FORMAT(' not possible to compute, first iteration',/)             40.41
!
!       *** number of points in which the penta-diagonal solver ***
!       *** did not converged                                   ***
!
        IF ( ITEST.GE.30 .OR. IDEBUG.EQ.1 ) THEN
          IF ((DYNDEP .OR. ICUR.EQ.1) .AND. INOCNV .NE. 0) THEN           40.31 40.23
             WRITE(PRINTF,122) INOCNV                                     40.23
          END IF
        END IF
 122    FORMAT(2X,'SIP solver: no convergence in ',I4,' gridpoints')      40.23
!
!----------------------------------------------------------------------   40.22
!     End master thread region.                                           40.22
!----------------------------------------------------------------------   40.22
!$OMP END MASTER                                                          40.22
!

!----------------------------------------------------------------------   40.22
!     Synchronize threads before checking accuracy and before             40.22
!     starting next iteration.                                            40.22
!----------------------------------------------------------------------   40.22
!$OMP BARRIER                                                             40.22
!
!       *** if accuracy has been reached then the iteration ***
!       *** can be terminated ---> goto 470                 ***
!
        IF ( (ITER.NE.1 .OR. PNUMS(21).EQ.0.) .AND.                       40.41
     &                                      ACCUR.GE.PNUMS(4) ) GOTO 470  40.41
!
 450  CONTINUE                                                            30.00
 470  CONTINUE
!TIMG      CALL SWTSTO(103)                                                    40.23

!----------------------------------------------------------------------   40.22
!     Begin deallocate private arrays.                                    40.22
!----------------------------------------------------------------------   40.22
!TIMG      CALL SWTSTA(101)                                                    40.23
      DEALLOCATE(IDCMIN)                                                  40.22
      DEALLOCATE(IDCMAX)                                                  40.22
      DEALLOCATE(ISCMIN)                                                  40.22
      DEALLOCATE(ISCMAX)                                                  40.22
      DEALLOCATE(CGO)                                                     40.22
      DEALLOCATE(KWAVE)                                                   40.22
      DEALLOCATE(CAX)                                                     40.22
      DEALLOCATE(CAY)                                                     40.22
      DEALLOCATE(CAS)                                                     40.22
      DEALLOCATE(CAD)                                                     40.22
      DEALLOCATE(CAX1)                                                    40.22
      DEALLOCATE(CAY1)                                                    40.22
      DEALLOCATE(ALIMW)                                                   40.22
      DEALLOCATE(UE)                                                      40.22
      DEALLOCATE(SA1)                                                     40.22
      DEALLOCATE(SA2)                                                     40.22
      DEALLOCATE(SFNL)                                                    40.22
      DEALLOCATE(DA1C)                                                    40.22
      DEALLOCATE(DA1P)                                                    40.22
      DEALLOCATE(DA1M)                                                    40.22
      DEALLOCATE(DA2C)                                                    40.22
      DEALLOCATE(DA2P)                                                    40.22
      DEALLOCATE(DA2M)                                                    40.22
      DEALLOCATE(DSNL)                                                    40.22
      DEALLOCATE(OBREDF)                                                  40.22
      DEALLOCATE(REFLSO)                                                  40.41
      DEALLOCATE(GROWW)                                                   40.22
      DEALLOCATE(ANYWND)                                                  40.22
      DEALLOCATE(SWMATR)                                                  40.22
      DEALLOCATE(LSWMAT)                                                  40.22
!TIMG      CALL SWTSTO(101)                                                    40.23
!----------------------------------------------------------------------   40.22
!     End deallocate private arrays.                                      40.22
!----------------------------------------------------------------------   40.22

!----------------------------------------------------------------------   40.22
!     End parallel region.                                                40.22
!----------------------------------------------------------------------   40.22
!$OMP END PARALLEL                                                        40.22
!
!     Print message when the solver did not converge in setup calculation
!
      IF (.NOT.CSETUP) THEN                                               30.82
          WRITE(PRINTF,124)                                               30.82
          IF (SCREEN.NE.PRINTF.AND.NSTATC.EQ.0.AND.IAMMASTER)             40.95 40.30 40.13
     &    WRITE(SCREEN,124)                                               30.82
      END IF                                                              30.82
 124  FORMAT(1X,'no convergence in set-up calculation')                   30.82
!
!TIMG      CALL SWTSTA(105)                                                    40.23
      IF (NPTST.GT.0 .AND. NSTATM.EQ.1                                    40.00
     &                    ) THEN
         IF (IFPAR.GT.0) WRITE (IFPAR, 11) CHTIME                         40.00
         IF (IFS1D.GT.0) WRITE (IFS1D, 11) CHTIME                         40.00
         IF (IFS2D.GT.0) WRITE (IFS2D, 11) CHTIME                         40.00
  11     FORMAT (A, T41, 'date-time')

         CALL PLTSRC (SWTSDA(1,1,1,JPWNDA)  ,SWTSDA(1,1,1,JPWNDB)  ,
     &                SWTSDA(1,1,1,JPWCAP)  ,SWTSDA(1,1,1,JPBTFR)  ,
     &                SWTSDA(1,1,1,JPWBRK)  ,SWTSDA(1,1,1,JP4S)    ,
     &                SWTSDA(1,1,1,JP4D)    ,SWTSDA(1,1,1,JPTRI)   ,
     &                SWTSDA(1,1,1,JPVEGT)  ,                             40.55
     &                AC2                   ,SPCSIG                ,      40.00
     &                COMPDA(1,JDP2)        ,XYTST                 ,
     &                                       KGRPNT                )      40.00
      END IF
!TIMG      CALL SWTSTO(105)                                                    40.23

!----------------------------------------------------------------------   40.22
!     Begin deallocate shared arrays.                                     40.22
!----------------------------------------------------------------------   40.22
!TIMG      CALL SWTSTA(101)                                                    40.23
      DEALLOCATE(SETPDA)                                                  40.22
      DEALLOCATE(HSAC1)                                                   40.22
      DEALLOCATE(HSAC2)                                                   40.22
      DEALLOCATE(SACC1)                                                   40.22
      DEALLOCATE(SACC2)                                                   40.22
      DEALLOCATE(HSAC0)                                                   40.41
      DEALLOCATE(HSDIFC)                                                  40.41
      DEALLOCATE(SACC0)                                                   40.93
      DEALLOCATE(TMDIFC)                                                  40.93
      DEALLOCATE(ISLMIN)                                                  40.23
      DEALLOCATE(NFLIM)                                                   40.23
      DEALLOCATE(NRSCAL)                                                  40.23
      DEALLOCATE(MEMNL4)                                                  40.22
!$    DEALLOCATE(LLOCK)                                                   40.22
      DEALLOCATE(AC2LOC)                                                  40.30
      DEALLOCATE(SWTSDA)                                                  40.31
!TIMG      CALL SWTSTO(101)                                                    40.23
!----------------------------------------------------------------------   40.22
!     End deallocate shared arrays.                                       40.22
!----------------------------------------------------------------------   40.22
!
      RETURN
      END subroutine SWCOMP
!
!************************************************************************
!
      SUBROUTINE SWOMPU (SWPDIR   ,KSX      ,KSY      ,
     &                   IX       ,IY       ,DDX      ,
     &                   DDY      ,DT       ,SNLC1    ,
     &                   DAL1     ,DAL2     ,DAL3     ,
     &                   XIS      ,SWTSDA   ,INOCNV   ,
     &                   AC2      ,COMPDA   ,SPCDIR   ,
     &                   SPCSIG   ,XYTST    ,ITER     ,                   30.72
     &                             CGO      ,                             40.41 40.22
     &                   CAX      ,CAY      ,CAS      ,
     &                   CAD      ,SWMATR   ,LSWMAT   ,                   30.90
     &                   KWAVE    ,
     &                   ALIMW    ,GROWW              ,                   40.17
     &                   UE       ,SA1      ,SA2      ,
     &                   DA1C     ,DA1P     ,DA1M     ,
     &                   DA2C     ,DA2P     ,DA2M     ,
     &                   SFNL     ,DSNL     ,MEMNL4   ,
     &                   IDCMIN   ,IDCMAX   ,
     &                   WWINT    ,WWAWG    ,WWSWG    ,
     &                   ISCMIN   ,ISCMAX   ,                             40.41
     &                   ANYWND   ,AC1      ,IT       ,
     &                   XCGRID   ,YCGRID   ,                             40.41 30.72
     &                   KGRPNT   ,CROSS    ,                             16/MAY
     &                   OBREDF   ,REFLSO   ,                             40.41 040697
     &                   ISLMIN   ,NFLIM    ,NRSCAL   ,                   40.23
     &                   CAX1,CAY1                                        33.08
     &                                                 )
!
!************************************************************************

      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
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
!  0. Authors
!
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.80: Nico Booij
!     30,81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     33.08: W. Erick Rogers (some S&L scheme-related changes)
!     33.09: Nico Booij (spherical coord.)
!     33.10: Nico Booij and Erick Rogers (2nd order upwind)
!     34.01: Jeroen Adema
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.08: Erick Rogers
!     40.09: Annette Kieftenburg
!     40.16: IJsbrand Haagsma
!     40.17: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!     40.23: Marcel Zijlema
!     40.28: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.61: Roop Lalbeharry
!
!  1. Updates
!
!     30.72, Nov. 97: Declaration of ISTAT, ITFRE, DDIR, DX, DY, GRAV,
!                     PI, U10 and WDIC removed because they are
!                     common and already declared in the INCLUDE file
!     32.02, Jan. 98: Introduced 1D-version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.72, Feb. 98: Modified argument list for update CGSTAB solver
!     30.70, Feb. 98: argument list of WINDP1 changed, current vel. added
!     40.00, July 98: KCGRD removed from Call WINDP1
!     40.00, Aug. 98: argument OBREDF added in call SPREDT
!                     subr SWTRCF called to calculate obstacle reduction factors
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description several variables
!     30.80, Nov. 98: Provision for limitation on Ctheta (refraction)
!     30.81, Jan. 99: Replaced variable STATUS by IERR (because STATUS is a
!                     reserved word)
!     34.01, Feb. 99: Introducing STPNOW
!     33.08  July 98: some S&L scheme-related changes
!     30.80, Aug. 99: Argument list SPROSD modified
!     40.10, Nov. 99: two arguments in call SWTRCF changed;
!                     CHS -> COMPDA(1,JHS) and WLEV2 -> COMPDA(1,JWLV2)
!     33.09, Nov. 99: call DSPHER added (ray curvature due to spherical coord.)
!     33.10, Jan. 00: changes re: the SORDUP scheme
!     40.09, May  00: Argument list SWTRCF modified
!     40.03, Jun. 00: new version of SPROSD has new argument list
!                     Ursell array added to argument list of SDISPA and SOURCE
!                     readability of test output improved
!     40.10, Sep. 00: Replaced SDISPA with SINTGRL
!     40.02, Sep. 00: Replaced SWMATR(1,1,JABIN) with LSWMAT (logical equivalence)
!     40.13, Mar. 01: if point was dry at previous time level, fall back to BSBT scheme
!                     order of calling SWAPAR and SPROXY changed, in order
!                     to get correct values of CGO as input to SPROXY
!     40.22, Sep. 01: Removed WAREA constructs and split SWMATR into      40.22
!                     SWMATR(real) and LSWMAT(logical).                   40.22
!     40.22, Sep. 01: Changed array definitions to use the parameter      40.22
!                     MICMAX instead of ICMAX.                            40.22
!     40.13, Oct. 01: loop over IC moved to subroutines SWAPAR and SPROXY
!     40.16, Dec. 01: Implementation of limiter switches
!     40.17, Dec. 01: Implementation of Multiple DIA
!     40.28, Dec. 01: Argument list SWTRCF modified                       40.28
!     40.23, Aug. 02: Print of CPU times added
!     40.23, Aug. 02: Introducing arrays NFLIM and NRSCAL
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.41, Aug. 04: code optimization
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Nov. 06: Hersbach and Janssen (1999) limiter option added
!
!  2. Purpose
!
!     This subroutine computes the wave spectrum for one sweep
!     direction, and is called four times per iteration.
!
!  3. Method
!
!
!    THIS IS THE STENCIL USED WITH THE S&L SCHEME:                  33.08
!                                                                   33.08
!      IY+1                      o 4                                33.08
!                                |                                  33.08
!                 8    6    2    | 1   5                            33.08
!      IY         O----O----O----*----O                             33.08
!                                |                                  33.08
!                           10   |                                  33.08
!      IY-1                 O----O 3                                33.08
!                                |                                  33.08
!                                |                                  33.08
!      IY-2                      O 7                                33.08
!                                |                                  33.08
!                                |                                  33.08
!      IY-3                      O 9                                33.08
!
!                 ^    ^    ^    ^    ^
!                 |    |    |    |    |                             33.08
!               IX-3 IX-2 IX-1  IX  IX+1
!
!     1: IX  , IY                                                   33.10
!     2: IX-1, IY                                                   33.10
!     3: IX  , IY-1                                                 33.10
!     4: IX  , IY+1                                                 33.10
!     5: IX+1, IY                                                   33.10
!     6: IX-2, IY                                                   33.10
!     7: IX  , IY-2                                                 33.10
!     8: IX-3, IY                                                   33.10
!     9: IX  , IY-3                                                 33.10
!    10: IX-1, IY-1                                                 33.10
!
!    THIS IS THE STENCIL USED WITH THE SORDUP SCHEME:               33.10
!                                                                   33.10
!                      4    2                                       33.10
!      IY              O----O----* 1                                33.10
!                                |                                  33.10
!                                |                                  33.10
!      IY-1                      O 3                                33.10
!                                |                                  33.10
!                                |                                  33.10
!      IY-2                      O 5                                33.10
!                                                                   33.10
!                      ^    ^    ^                                  33.10
!                      |    |    |                                  33.10
!                    IX-2 IX-1  IX                                  33.10
!     1: IX  , IY                                                   33.10
!     2: IX-1, IY                                                   33.10
!     3: IX  , IY-1                                                 33.10
!     4: IX-2, IY                                                   33.10
!     5: IX  , IY-2                                                 33.10
!
!  4. Argument variables
!
!     ITER  : input Iteration counter for SWAN
!     IT    : input Time step counter for SWAN
!
      INTEGER ITER,   IT
!
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
      REAL :: SPCDIR(MDC,6)                                               30.82
      REAL :: SPCSIG(MSC)                                                 30.72
      REAL :: XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!     Since the real piece of LSWMAT was removed, the dimensions of       40.22
!     LSWMAT are now MLSWMAT not MSWMAT.                                  40.22
      LOGICAL :: LSWMAT(MDC,MSC,MLSWMAT)                                  30.90 40.22
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     SPREDT
!     SWAPAR
!     SPROXY
!     SINTGRL
!     SOURCE
!     ACTION
!     SOLMAT
!     SOLMT1
!     SPROSD
!     SWPSEL
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     ---------------------------------------------------------
!     Call WINDP2 to compute some wave parameters necessarry for
!                 the wind subroutines. The wind sea energy spectrum
!                 is computed before every iteration
!     Compute for the two nearby points:
!       {to reduce the size of the arrays K, CPX, CPY, CAX, CAY, CAS, CAD
!       and CGO, CP use a FUNCTION ICODE(_,_) in were the information
!       of the nearby gridpoints is stored.
!       The size of the arrays of the wave parameters are reduced significantly,
!       par example: CAX(ID,IS,IX,IY) --> CAX(ID,IS,ICMAX)  with ICMAX = 3
!       If a higher order scheme is used ICMAX can be increased so that points
!       at locations ksx = -2,+2 and ksy = -2,+2 can be used:
!
!                                 o ksy=+2
!                                 |
!                                 |
!                                 o ksy=+1     with:  * = (0,0)
!                                 |
!                                 |
!                   o------o------*------o------o
!                ksx=-2   ksx=-1  |   ksx=+1   ksx=+2
!                                 |
!                                 o ksy=-1
!                                 |
!                                 |
!                                 o ksy=-2
!
!          Molecule:
!
!          (4)     2
!            o------o------* 1           Central grid point     : IC = 1              30.70(?), 33.10
!                          |             Point in X-direction   : IC = 2              30.70(?), 33.10
!                          |             Point in Y-direction   : IC = 3              30.70(?), 33.10
!                        3 o             Point in X-direction   : IC = (4)            30.70(?), 33.10
!                          |             Point in Y-diretion    : IC = (5)            30.70(?), 33.10
!                          |             5 gridpoints --> ICC = 5                     30.70(?), 33.10
!                      (5) o             ( ) = is not used by default (BSBT) scheme   30.70(?), 33.10
!
!          Notice that IX and IY are still in the argument list because of
!          the counter of DEP2(IX,IY) and UX2(IX,IY) and UY2(IX,IY) !
!
!     For every ICC = 1 to ICMAX do
!       If (ICC = 1) then
!         IC  = ICODE(0,0)                {central gridpoint}
!         ICX = IX
!         ICY = IY
!         --------------------------------------------------
!       Else if (ICC = 2) then
!         IC = ICODE(KSX,0)               {left or right gridpoint}
!         ICX = IX+KSX
!         ICY = IY
!         --------------------------------------------------
!       Else if (ICC = 3 ) then
!         IC = ICODE(0,KSY)               {top or bottom gridpoint)
!         ICX = IX
!         ICY = IY+KSY
!         --------------------------------------------------
!     End if
!
!     -----------------------------------------------------------------
!     For each gridpoint (IC=1,2,3)  do:
!       Call SWAPAR  to compute the wavenumber K and the group velocity CGO
!       Call SPROXY  to compute propagation velocities CAX, CAY of
!                    energy propagates
!       If central gridpoint (IC=1) do
!         Call SWGEOM  to compute geometric quantities due to
!                      curvilinear grid
!         Call SWPSEL  to compute the bins that fall within a sweep
!                      and which are propagated within a sweep
!         Call SPROSD  to calculate the propagation velocities in
!                      spectral space (CAS and CAD)
!     -------------------------------------------------------
!     If depth > 1.e-4 then do
!       If wind is present :
!         Call WINDP1 to compute the wind speed, PM frequency, mean
!                     wind direction, wind friction velocity and counters
!       -----------------------------------------------------------------
!       Call CNTAIL to compute the contribution of high frequency
!                   tail to the spectrum
!       -----------------------------------------------------------------
!       For the first iteration (ITER = 1), do:
!         Call SPREDT to estimate the action density in stationary mode
!       -----------------------------------------------------------------
!       Call SINTGRL to compute some wave parameters (mean frequency
!                    mean wave number, near bottom velocity and signi-
!                    ficant wave height, fraction of breaking waves
!       -----------------------------------------------------------------
!       Call SOURCE  to compute the source terms for each bin which fall
!                    within a sweep:
!                   1. Dissipation by wave-bottom effects
!                   2. Dissipation due to surf breaking
!                   3. Dissipation due to whitecapping
!                   4. Generation of wave energy by wind effects
!                   5. Nonlinear wave-wave interactions (quadruplets)
!                   6. Nonlinear wave-wave interactions (triads)
!                   7. Dissipation due to vegetation
!       -----------------------------------------------------------------
!       Call ACTION  calculate the derivatives in x,y,s,d space and store
!                    the results in the corresponding arrays
!       -----------------------------------------------------------------
!       If a current is present do
!         If implicit scheme in frequency space do
!         -----
!           Call SWSIP    to solve penta-diagonal system by means
!                         of Stone's SIP solver
!         -----
!         endif
!       else if no current is present do
!         Call SOLMAT   to solve tri-diagonal system
!       end if
!     ------------------------------------------------------------------
!     CALL WINDP3   for a first or second generation model: limit the
!                   computed action density in a gridpoint according to
!                   the saturation spectra.
!     ---------------------------------------------------------
!     End of SWOMPU
!     ---------------------------------------------------------
!
! 13. Source text
!
      INTEGER  IC    ,IX    ,IY    ,IS    ,SWPDIR,                        40.00
     &         KSX   ,KSY   ,                                             40.00
     &         IDWMIN,IDWMAX,IDTOT ,ISTOT ,IDDLOW,IDDTOP,
     &         ISSTOP,INOCNV                                              40.00
      INTEGER  LINK(MICMAX)                                               33.09
!
      REAL     DDX   ,DDY   ,DT    ,                                      40.00
     &         ETOT  ,AC2TOT,ABRBOT,HM    ,HS    ,QBLOC ,                 40.00
     &         SMESPC,KMESPC,ETOTW ,WIND10,FPM   ,
     &         THETAW,SNLC1 ,DAL1  ,DAL2  ,DAL3  ,XIS   ,
     &         UFRIC ,SMEBRK
!
      LOGICAL  INSIDE                                                     33.09
      LOGICAL  LPREDT
!
      INTEGER :: XYTST(2*NPTST) ,IDCMIN(MSC)                              40.22
      INTEGER :: IDCMAX(MSC)    ,ISCMIN(MDC)    ,ISCMAX(MDC)
      INTEGER :: WWINT(*)
      INTEGER :: KGRPNT(MXC,MYC)                                          40.00
      INTEGER :: CROSS(2,MCGRD)
!
!     *** number of arrays for SWAN ***
!
      REAL  :: AC2(MDC,MSC,MCGRD)                                         30.21
      REAL  :: AC1(MDC,MSC,MCGRD)                                         30.00
      REAL  :: COMPDA(MCGRD,MCMVAR)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  :: CGO(MSC,MICMAX)          ,                                 40.22
     &         CAX(MDC,MSC,MICMAX)      ,                                 40.22
     &         CAY(MDC,MSC,MICMAX)      ,                                 40.22
     &         CAX1(MDC,MSC,MICMAX)     ,                                 33.08 40.22
     &         CAY1(MDC,MSC,MICMAX)     ,                                 33.08 40.22
     &         CAS(MDC,MSC,MICMAX)      ,                                 40.22
     &         CAD(MDC,MSC,MICMAX)                                        40.22
      REAL  :: ALIMW(MDC,MSC)
!              Since the logical piece of SWMATR was removed, the         40.22
!              dimensions of SWMATR are now MSWMATR not MSWMAT            40.22
      REAL  :: SWMATR(MDC,MSC,MSWMATR)                                    40.22
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  :: KWAVE(MSC,MICMAX),                                         40.17 40.22
     &         UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,
     &         SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &         SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &         DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &         MEMNL4(MDC,MSC,MCGRD)               ,                      30.21
     &         SWTSDA(MDC,MSC,NPTSTA,MTSVAR)         ,                    40.00
     &         WWAWG(*)                              ,
     &         WWSWG(*)                              ,
     &         RDX(10)       ,RDY(10)               ,                     15/MAY 40.08
     &         OBREDF(MDC,MSC,2)                    ,                     040697
     &         REFLSO(MDC,MSC)                                            40.41
!
      LOGICAL  GROWW(MDC,MSC)    ,
     &         ANYWND(MDC)
!
      INTEGER :: ISLMIN(MCGRD), NFLIM(MCGRD), NRSCAL(MCGRD)               40.23
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWOMPU')
!     *** Get grid point numbers for points in computational stencil ***
      IXCGRD(1) = IX                                                      40.00
      IYCGRD(1) = IY                                                      40.00
      KCGRD(1)  = KGRPNT(IX,IY)                                           30.21
      IF (KCGRD(1).GT.1) THEN                                             40.00
        IXCGRD(2) = IX+KSX                                                40.00
        IYCGRD(2) = IY                                                    40.00
        IXCGRD(3) = IX                                                    40.00
        IYCGRD(3) = IY+KSY                                                40.00
        PROPSL = PROPSC                                                   33.09
        IF (PROPSC.EQ.1) THEN                                             33.09
          ICMAX = 3
        ELSE IF (PROPSC.EQ.2) THEN          ! SORDUP scheme               33.10
!         add more points for higher order scheme                         33.10
          ICMAX = 5                                                       33.10
          IXCGRD(4) = IX+2*KSX                                            33.10
          IYCGRD(4) = IY                                                  33.10
          IXCGRD(5) = IX                                                  33.10
          IYCGRD(5) = IY+2*KSY                                            33.10
        ELSE IF (PROPSC.EQ.3) THEN          ! S&L scheme                  33.09
!         add more points for higher order schemes here
          ICMAX = 10                                                      33.09
          IXCGRD(4) = IX                                                  33.09
          IYCGRD(4) = IY-KSY                                              33.09
          IXCGRD(5) = IX-KSX                                              33.09
          IYCGRD(5) = IY                                                  33.09
          IXCGRD(6) = IX+2*KSX                                            33.09
          IYCGRD(6) = IY                                                  33.09
          IXCGRD(7) = IX                                                  33.09
          IYCGRD(7) = IY+2*KSY                                            33.09
          IXCGRD(8) = IX+3*KSX                                            33.09
          IYCGRD(8) = IY                                                  33.09
          IXCGRD(9) = IX                                                  33.09
          IYCGRD(9) = IY+3*KSY                                            33.09
          IXCGRD(10) = IX+KSX                                             33.09
          IYCGRD(10) = IY+KSY                                             33.09
        ENDIF
        DO IC = 2, ICMAX                                                  40.00
!         if one of the points of a stencil is outside the computational  33.09
!         domain, fall back to first order scheme                         33.09
          INSIDE = .TRUE.                                                 33.09
          IF (IXCGRD(IC).LT.1) THEN                                       33.09
            IF (KREPTX.GT.0) THEN                                         33.09
!             domain is repeating in x-direction                          33.09
              IXCGRD(IC) = IXCGRD(IC) + MXC                               33.09
            ELSE                                                          33.09
              IF (IC.LE.3) THEN                                           33.09
                PROPSL = 0                                                33.09
              ELSE                                                        33.09
                IF (PROPSC.GT.1) PROPSL = 1                               33.09
              ENDIF                                                       33.09
              INSIDE = .FALSE.                                            33.09
            ENDIF                                                         33.09
          ENDIF                                                           33.09
          IF (IXCGRD(IC).GT.MXC) THEN                                     33.09
            IF (KREPTX.GT.0) THEN                                         33.09
              IXCGRD(IC) = IXCGRD(IC) - MXC                               33.09
            ELSE                                                          33.09
              IF (IC.LE.3) THEN                                           33.09
                PROPSL = 0                                                33.09
              ELSE                                                        33.09
                IF (PROPSC.GT.1) PROPSL = 1                               33.09
              ENDIF                                                       33.09
              INSIDE = .FALSE.                                            33.09
            ENDIF                                                         33.09
          ENDIF                                                           33.09
          IF (IYCGRD(IC).LT.1 .OR. IYCGRD(IC).GT.MYC) THEN                33.09
            IF (.NOT.ONED) THEN                                           33.09
              IF (IC.LE.3) THEN                                           33.09
                PROPSL = 0                                                33.09
              ELSE                                                        33.09
                IF (PROPSC.GT.1) PROPSL = 1                               33.09
              ENDIF                                                       33.09
            ENDIF                                                         33.09
            INSIDE = .FALSE.                                              33.09
          ENDIF                                                           33.09
          IF (INSIDE) THEN                                                33.09
            KCGRD(IC) = KGRPNT(IXCGRD(IC),IYCGRD(IC))                     40.00
          ELSE                                                            33.09
            KCGRD(IC) = 1                                                 33.09
          ENDIF                                                           33.09
!         if point in stencil is dry, fall back to BSBT scheme.           33.10
          IF (PROPSC.GT.1) THEN                                           33.10
            IF (COMPDA(KCGRD(IC),JDP2).LE.DEPMIN) PROPSL = 1              33.10
            IF (NSTATC.GT.0) THEN                                         40.13
!             if nonstationary, check also previous time level            40.13
              IF (COMPDA(KCGRD(IC),JDP1).LE.DEPMIN) PROPSL = 1            40.13
            ENDIF                                                         40.13
          END IF                                                          33.10
        ENDDO
        IF (PROPSL.EQ.0) THEN                                             33.09
          ICMAX = 1                                                       33.09
        ELSEIF (PROPSL.EQ.1) THEN                                         33.09
          ICMAX = 3                                                       33.09
        ENDIF                                                             33.09
      ELSE
        PROPSL = 0                                                        40.13
        ICMAX = 1                                                         40.13
      ENDIF
!
!     *** If there are obstacles crossing the points in the stencil ***
!     *** then fall back to first order scheme                      ***
!
!TIMG      CALL SWTSTA(136)                                                    40.23
      IF (NUMOBS.NE.0 .AND. PROPSL.NE.1) THEN
         IF (PROPSL.EQ.3) THEN
            NLINK  = 10
            IF (SWPDIR .EQ. 1 ) THEN
               LINK(1)  = CROSS(1,KCGRD(1))                               33.08
               LINK(2)  = CROSS(2,KCGRD(1))
               LINK(3)  = CROSS(1,KCGRD(2))
               LINK(4)  = CROSS(2,KCGRD(2))
               LINK(5)  = CROSS(1,KCGRD(3))
               LINK(6)  = CROSS(2,KCGRD(3))
               LINK(7)  = CROSS(2,KCGRD(4))
               LINK(8)  = CROSS(1,KCGRD(5))
               LINK(9)  = CROSS(1,KCGRD(6))
               LINK(10) = CROSS(2,KCGRD(7))
            ELSE IF (SWPDIR .EQ. 2) THEN
               LINK(1)  = CROSS(1,KCGRD(2))                               33.08
               LINK(2)  = CROSS(2,KCGRD(1))
               LINK(3)  = CROSS(1,KCGRD(6))
               LINK(4)  = CROSS(2,KCGRD(2))
               LINK(5)  = CROSS(1,KCGRD(10))
               LINK(6)  = CROSS(2,KCGRD(3))
               LINK(7)  = CROSS(2,KCGRD(4))
               LINK(8)  = CROSS(1,KCGRD(1))
               LINK(9)  = CROSS(1,KCGRD(8))
               LINK(10) = CROSS(2,KCGRD(7))
            ELSE IF (SWPDIR .EQ. 3) THEN
               LINK(1)  = CROSS(1,KCGRD(2))                               33.08
               LINK(2)  = CROSS(2,KCGRD(3))
               LINK(3)  = CROSS(1,KCGRD(6))
               LINK(4)  = CROSS(2,KCGRD(10))
               LINK(5)  = CROSS(1,KCGRD(10))
               LINK(6)  = CROSS(2,KCGRD(7))
               LINK(7)  = CROSS(2,KCGRD(1))
               LINK(8)  = CROSS(1,KCGRD(1))
               LINK(9)  = CROSS(1,KCGRD(8))
               LINK(10) = CROSS(2,KCGRD(9))
            ELSE IF (SWPDIR .EQ. 4) THEN
               LINK(1)  = CROSS(1,KCGRD(1))                               33.08
               LINK(2)  = CROSS(2,KCGRD(3))
               LINK(3)  = CROSS(1,KCGRD(2))
               LINK(4)  = CROSS(2,KCGRD(10))
               LINK(5)  = CROSS(1,KCGRD(3))
               LINK(6)  = CROSS(2,KCGRD(7))
               LINK(7)  = CROSS(2,KCGRD(1))
               LINK(8)  = CROSS(1,KCGRD(5))
               LINK(9)  = CROSS(1,KCGRD(6))
               LINK(10) = CROSS(2,KCGRD(9))
            ENDIF
         ELSE IF (PROPSL.EQ.2) THEN                                       33.10
            NLINK  = 4                                                    33.10
            IF (SWPDIR .EQ. 1 ) THEN                                      33.10
               LINK(1)  = CROSS(1,KCGRD(1))                               33.10
               LINK(2)  = CROSS(2,KCGRD(1))                               33.10
               LINK(3)  = CROSS(1,KCGRD(2))                               33.10
               LINK(4)  = CROSS(2,KCGRD(3))                               33.10
            ELSE IF (SWPDIR .EQ. 2) THEN                                  33.10
               LINK(1)  = CROSS(1,KCGRD(2))                               33.10
               LINK(2)  = CROSS(2,KCGRD(1))                               33.10
               LINK(3)  = CROSS(1,KCGRD(4))                               33.10
               LINK(4)  = CROSS(2,KCGRD(3))                               33.10
            ELSE IF (SWPDIR .EQ. 3) THEN                                  33.10
               LINK(1)  = CROSS(1,KCGRD(2))                               33.10
               LINK(2)  = CROSS(2,KCGRD(3))                               33.10
               LINK(3)  = CROSS(1,KCGRD(4))                               33.10
               LINK(4)  = CROSS(2,KCGRD(5))                               33.10
            ELSE IF (SWPDIR .EQ. 4) THEN                                  33.10
               LINK(1)  = CROSS(1,KCGRD(1))                               33.10
               LINK(2)  = CROSS(2,KCGRD(3))                               33.10
               LINK(3)  = CROSS(1,KCGRD(2))                               33.10
               LINK(4)  = CROSS(2,KCGRD(5))                               33.10
            ENDIF                                                         33.10
         ENDIF                                                            33.10
         IF (PROPSL.GT.1) THEN                                            33.10
!           if there is an obstacle crossing, fall back to 1st order      33.10
            DO ILINK = 1, NLINK                                           33.08
               IF (LINK(ILINK).GT.0) PROPSL = 1                           33.08
            ENDDO
         ENDIF
      ENDIF
!TIMG      CALL SWTSTO(136)                                                    40.23

      IF (ITEST .GE. 180 ) THEN
        WRITE(PRINTF,188) SWPDIR
 188    FORMAT(' Points in stencil in subr SWOMPU, sweep : ',I1,
     &  /,'POINT( IX, IY),  INDEX,    COORDX,       COORDY')
!
        DO IC = 1, ICMAX                                                  40.00
          WRITE(PRINTF,187) IXCGRD(IC), IYCGRD(IC), KCGRD(IC),            40.00
     &    XCGRID(IXCGRD(IC),IYCGRD(IC)), YCGRID(IXCGRD(IC),IYCGRD(IC))    40.00
 187      FORMAT(4X,I4,1X,I4,3X,I5,5X,F10.2,4X,F10.2)
        ENDDO                                                             40.00
      ENDIF
!
!     *** determine whether the point is a test point ***
!
      IPTST  = 0
      TESTFL = .FALSE.
      IF (NPTST.GT.0) THEN
        DO 20 II = 1, NPTST
          IF (IX.NE.XYTST(2*II-1)) GOTO 20
          IF (IY.NE.XYTST(2*II  )) GOTO 20
!
          IPTST = II
          TESTFL = .TRUE.
          IF (ITEST .GE. 10)
     &    WRITE(PRINTF, 18) IPTST, IX+MXF-2, IY+MYF-2, KCGRD(1), ITER,    40.30
     &                      SWPDIR                                        40.30
  18      FORMAT(' Test point ', I2, ', (ix,iy)', 2I5,
     &           ', point index ',I5,                                     30.21
     &           ', iter ', I2, ', sweep ', I1)
  20    CONTINUE
      ENDIF
!
      IF (TESTFL .AND. ITEST .GE. 220) THEN                               40.00
        INDEX = KGRPNT(IX,IY)
        WRITE(PRINTF,322) IX+MXF-1,IY+MYF-1,INDEX                         40.30
 322    FORMAT(' Action densities for IX  IY IDNX: ', 2I4, I7)            40.00
        DO IS = 1, MSC
          WRITE(PRINTF, 323) IS
 323      FORMAT(' frequency ', I4)
          WRITE(PRINTF, 324) (AC2(ID,IS,INDEX), ID=1, MDC)                40.00
 324      FORMAT (100E12.4)
        ENDDO
      ENDIF
!
      IF(TESTFL.AND.ITEST.GE.10) WRITE(PRINTF,2321) SWPDIR,IX+MXF-1,      40.30
     &                                              IY+MYF-1              40.30
 2321 FORMAT(//,' sweep direction and node IX, IY ',3I4)
!
!     --- in case of non-active point set action density in land points
!         equal to zero and go back to main program
      IF (KCGRD(1).LE.1 .OR. COMPDA(KCGRD(1),JDP2).LE.DEPMIN) THEN        40.41 25/MAR
         DO IS = 1, MSC
            DO ID = 1, MDC
               AC2(ID,IS,KCGRD(1)) = 0.                                   30.21
            ENDDO
         ENDDO
         COMPDA(KCGRD(1),JHS) = 0.                                        40.41
         RETURN
      END IF
!
      IF (KCGRD(2) .LE. 1) RETURN                                         30.70
      IF (.NOT.ONED .AND. (KCGRD(3).LE.1)) RETURN                         30.70
!
      IF (PROPSL.EQ.3 .AND. NSTATC.GT.0) THEN                             33.08

!         calculate propagation velocities for old time level             40.13
!         (needed for S&L scheme nonstationary)                           40.13

!         COMPDA(1,JDP2) is dep2 ;change to..dep1 which is COMPDA(1,JDP1) 33.08
!         COMPDA(1,JVX2) is ux2  ;change to..ux1 which is COMPDA(1,JVX1)  33.08
!         COMPDA(1,JVY2) is uy2  ;change to..uy1 which is COMPDA(1,JVY1)  33.08
!
!         we could save CPU time by calculating the CAX values only once  33.08
!         when CAX is constant, but this would require SWAN to            33.08
!         save CAX over the entire grid (more memory).                    33.08
!
!         if nonstationary, then CAX1 is not necessarily equal to CAX     33.08
!         also, if we are using the BSBT scheme only,                     33.08
!         then CAX1, CAY1 are not needed.                                 33.08
!
!TIMG          CALL SWTSTA(110)                                                40.23
          CALL SWAPAR ( COMPDA(1,JDP1), KWAVE, CGO, SPCSIG )              40.41
!TIMG          CALL SWTSTO(110)                                                40.23
!
!         *** compute the propagation velocities CAX1 and CAY1       ***
!         *** for all directions for the gridpoints IC = 1 to ICMAX  ***

!TIMG          CALL SWTSTA(111)                                                40.23
          CALL SPROXY (CAX1            ,
     &                 CAY1           ,CGO            ,SPCDIR(1,2)    ,   33.08
     &                 SPCDIR(1,3)    ,COMPDA(1,JVX1) ,COMPDA(1,JVY1) ,   33.08
     &                 SWPDIR
     &                                                              )     33.08
!TIMG          CALL SWTSTO(111)                                                40.23
!
      END IF                                                              33.08

!     calculate propagation velocities for new time level                 40.13

!     *** Compute wavenumber KWAVE and group velocity CGO   ***
!     *** in the gridpoints of the stencil                  ***

!TIMG      CALL SWTSTA(110)                                                    40.23
      CALL SWAPAR ( COMPDA(1,JDP2), KWAVE, CGO, SPCSIG )                  40.41
!TIMG      CALL SWTSTO(110)                                                    40.23
!
!     *** compute the propagation velocities CAX and CAY        ***
!     *** for all directions for the gridpoints IC = 1 to ICMAX ***
!
!TIMG      CALL SWTSTA(111)                                                    40.23
      CALL SPROXY (CAX            ,
     &             CAY            ,CGO            ,SPCDIR(1,2)    ,
     &             SPCDIR(1,3)    ,COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,
     &             SWPDIR
     &                                                            )
!TIMG      CALL SWTSTO(111)                                                    40.23
!
!     --- compute geometric quantities due to curvilinear grid            40.41
!
!TIMG      CALL SWTSTA(112)
      CALL SWGEOM ( RDX, RDY, XCGRID, YCGRID, SWPDIR )                    40.41
!TIMG      CALL SWTSTO(112)
!
!     *** compute minimum and maximum counter (IDCMIN and ***
!     *** IDCMAX) and fill the array ANYBIN to determine  ***
!     *** if a bin lies within the sweep considered       ***
!
!TIMG      CALL SWTSTA(112)                                                    40.23
      CALL SWPSEL (SWPDIR                            ,IDCMIN        ,
     &             IDCMAX          ,CAX              ,
     &             CAY             ,LSWMAT(1,1,JABIN),                    40.02
     &                              ISCMIN           ,                    40.00
     &             ISCMAX          ,IDTOT            ,ISTOT         ,
     &             IDDLOW          ,IDDTOP           ,ISSTOP        ,
     &             COMPDA(1,JDP2)  ,COMPDA(1,JVX2)   ,COMPDA(1,JVY2),
     &             SPCDIR          ,RDX              ,RDY           ,     40.41 30.72
     &             KGRPNT                                                 40.13
     &                                                              )
!TIMG      CALL SWTSTO(112)                                                    40.23

!     *** compute the propagation velocities CAS and CAD   ***
!     *** for the central gridpoint only and only for the  ***
!     *** directional domain : IDCMIN-1 until IDCMAX+1     ***

!TIMG      CALL SWTSTA(113)                                                    40.23
      CALL SPROSD (SPCSIG         ,KWAVE          ,CAS            ,       40.03
     &             CAD            ,CGO            ,                       30.80
     &             COMPDA(1,JDP2) ,COMPDA(1,JDP1) ,SPCDIR(1,2)    ,
     &             SPCDIR(1,3)    ,COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,
     &             SWPDIR         ,IDCMIN         ,IDCMAX         ,
     &             SPCDIR(1,4)    ,SPCDIR(1,6)    ,SPCDIR(1,5)    ,       30.80
     &             RDX            ,RDY            ,                       30.80
     &             CAX            ,CAY            ,LSWMAT(1,1,JABIN),     40.02
     &             KGRPNT         ,XCGRID         ,YCGRID         ,       40.03
     &             IDDLOW         ,IDDTOP                                 40.61
     &                                                            )
!TIMG      CALL SWTSTO(113)                                                    40.23

!TIMG      CALL SWTSTA(114)                                                    40.23
      IF (KSPHER.GT.0 .AND. IREFR.NE.0) THEN                              40.41
!
!        *** compute the change of propagation velocity CAD   ***
!        *** due to the use of spherical coordinates          ***
!
         CALL DSPHER (CAD                ,CAX            ,                40.41 33.09
     &                CAY                ,                                40.41
     &                LSWMAT(1,1,JABIN)  ,YCGRID         ,                40.02
     &                SPCDIR(1,2)        ,SPCDIR(1,3)    )                40.41
      ENDIF
!TIMG      CALL SWTSTO(114)                                                    40.23
!
      IF ( IDTOT.GT.0 ) THEN                                              40.41 25/MAR
!
!       *** initialize friction velocity and Fpm frequency even ***
!       *** when there is no wind input                         ***
!
        UFRIC = 1.E-15
        FPM   = 1.E-15
        IF ( IWIND .GE. 1 ) THEN
!
!         *** compute the wind speed, mean wind direction, the    ***
!         *** PM frequency, wind friction velocity U*  and the    ***
!         *** minimum and maximum counters for active wind input  ***
!
!TIMG          CALL SWTSTA(115)                                                40.23
          CALL WINDP1 (WIND10     ,THETAW     ,
     &                 IDWMIN     ,IDWMAX     ,
     &                 FPM        ,UFRIC      ,
     &                 COMPDA(1,JWX2) ,COMPDA(1,JWY2) ,
     &                 ANYWND     ,SPCDIR     ,                           40.00
     &                 COMPDA(1,JVX2) ,COMPDA(1,JVY2) ,SPCSIG     )       30.70
!TIMG          CALL SWTSTO(115)                                                40.23
        END IF
!
!       *** estimate action density in case of first iteration ***
!       *** in stationary mode (since it is zero in first      ***
!           stationary run)                                    ***
!
        LPREDT=.FALSE.
        IF (ICOND.NE.4 .AND. ITER.EQ.1 .AND. NSTATC.EQ.0) THEN            40.41 40.00
          LPREDT=.TRUE.
          COMPDA(KCGRD(1),JHS) = 0.
          GOTO 333
        END IF
 330    IF (LPREDT) THEN
          CALL SPREDT (SWPDIR           ,AC2               ,CAX       ,
     &                 CAY              ,IDCMIN            ,IDCMAX    ,
     &                 ISSTOP           ,LSWMAT(1,1,JABIN) ,              40.02
     &                 RDX              ,RDY               ,OBREDF    )   40.00
           LPREDT=.FALSE.
        END IF
!
!       Calculate various integral parameters for use in the source terms
!
!TIMG        CALL SWTSTA(116)                                                  40.23
        CALL SINTGRL  (SPCDIR  ,KWAVE   ,AC2     ,                        40.02
     &                 COMPDA(1,JDP2)   ,QBLOC   ,COMPDA(1,JURSEL),       40.02
     &                 RDX     ,RDY     ,                                 40.02
     &                 AC2TOT  ,ETOT    ,                                 40.02
     &                 ABRBOT  ,COMPDA(1,JUBOT)  ,HS      ,               40.02
     &                 COMPDA(1,JQB)    ,                                 40.02
     &                 HM      ,KMESPC  ,SMEBRK,                          40.02
     &                 COMPDA(1,JPBOT)  ,                                 40.51
     &                 SWPDIR                  )                          40.16
!TIMG        CALL SWTSTO(116)                                                  40.23
!
        COMPDA(KCGRD(1),JHS) = HS                                         30.70
!
!       *** If there are obstacles crossing the points in the stencil ***
!       *** then the transmission and reflection coeff. are computed  ***
!       *** and also the contribution to the source term              ***
!
 333    CONTINUE
!TIMG        CALL SWTSTA(136)                                                  40.23
        IF (NUMOBS .NE. 0) THEN
!
!         *** OBREDF(:,:,2) are the transmission coeff for the two links ***
!         *** in the stencil (between the three point on the stencil)    ***
!         *** REFLSO(:,:) contains the contribution to the source term   ***
!
          DO IDC = 1, MDC
            DO ISC = 1, MSC
              OBREDF(IDC,ISC,1) = 1.                                      40.00
              OBREDF(IDC,ISC,2) = 1.                                      40.00
              REFLSO(IDC,ISC  ) = 0.
            ENDDO
          ENDDO
!
          IF (SWPDIR .EQ. 1 ) THEN
             LINK(1) = CROSS(1,KCGRD(1))
             LINK(2) = CROSS(2,KCGRD(1))
          ELSE IF (SWPDIR .EQ. 2) THEN
             LINK(1) = CROSS(1,KCGRD(2))
             LINK(2) = CROSS(2,KCGRD(1))
          ELSE IF (SWPDIR .EQ. 3) THEN
             LINK(1) = CROSS(1,KCGRD(2))
             LINK(2) = CROSS(2,KCGRD(3))
          ELSE IF (SWPDIR .EQ. 4) THEN
             LINK(1) = CROSS(1,KCGRD(1))
             LINK(2) = CROSS(2,KCGRD(3))
          ENDIF

          IF (LINK(1) .NE. 0 .OR. LINK(2) .NE. 0) THEN                    40.00
             IF (ITEST .GE. 120) WRITE(PRINTF,10)
     &          SWPDIR,KCGRD(1),LINK(1),LINK(2)
 10          FORMAT(' SWOMPU:  SWPDIR POINT LINK1  LINK2 = ',4(1X,I5))
!
             CALL SWTRCF (COMPDA(1,JWLV2),                                40.31 40.09
     &                    COMPDA(1,JHS), LINK, OBREDF,                    40.09
     &                    AC2, REFLSO, KGRPNT, XCGRID,                    40.41 40.09
     &                    YCGRID, CAX, CAY, RDX, RDY, LSWMAT(1,1,JABIN),  40.02
     &                    SPCSIG, SPCDIR)                                 40.13 40.28
          ENDIF
!
        ENDIF
!TIMG        CALL SWTSTO(136)                                                  40.23
        IF (LPREDT) GOTO 330
!
!       *** compute source terms and fill the matrix ***
!
!TIMG        CALL SWTSTA(117)                                                  40.23
        CALL SOURCE (ITER   ,IX                  ,IY                  ,
     &  SWPDIR              ,KWAVE               ,SPCSIG              ,   30.72
     &  SPCDIR(1,2)         ,SPCDIR(1,3)         ,AC2                 ,
     &  COMPDA(1,JDP2)      ,SWMATR(1,1,JMATD)   ,SWMATR(1,1,JMATR)   ,
     &  ABRBOT              ,KMESPC              ,SMESPC              ,
     &  COMPDA(1,JUBOT)     ,UFRIC               ,COMPDA(1,JVX2)      ,
     &  COMPDA(1,JVY2)      ,IDCMIN              ,IDCMAX              ,
     &  IDDLOW              ,IDDTOP              ,IDWMIN              ,
     &  IDWMAX              ,ISSTOP              ,SWTSDA(1,1,1,JPWNDA),
     &  SWTSDA(1,1,1,JPWNDB),SWTSDA(1,1,1,JPWCAP),SWTSDA(1,1,1,JPBTFR),
     &  SWTSDA(1,1,1,JPWBRK),SWTSDA(1,1,1,JP4S)  ,SWTSDA(1,1,1,JP4D)  ,
     &  SWTSDA(1,1,1,JPVEGT),                                             40.55
     &  SWTSDA(1,1,1,JPTRI) ,                     HS                  ,   40.22
     &  ETOT                ,QBLOC               ,THETAW              ,
     &  HM                  ,FPM                 ,WIND10              ,
     &  ETOTW               ,GROWW               ,ALIMW               ,
     &  SMEBRK              ,SNLC1               ,
     &  DAL1                ,DAL2                ,DAL3                ,
     &  UE                  ,SA1                 ,                        40.17
     &  SA2                 ,DA1C                ,DA1P                ,
     &  DA1M                ,DA2C                ,DA2P                ,
     &  DA2M                ,SFNL                ,DSNL                ,
     &  MEMNL4              ,WWINT               ,WWAWG               ,
     &  WWSWG               ,CGO                 ,COMPDA(1,JUSTAR)    ,
     &  COMPDA(1,JZEL)      ,SPCDIR              ,ANYWND              ,
     &  SWMATR(1,1,JDIS0)   ,SWMATR(1,1,JDIS1)   ,
     &  SWMATR(1,1,JGEN0)   ,SWMATR(1,1,JGEN1)   ,                        40.85
     &  SWMATR(1,1,JRED0)   ,SWMATR(1,1,JRED1)   ,                        40.85
     &  XIS                 ,COMPDA(1,JFRC2)     ,IT                  ,   40.00
     &  COMPDA(1,JNPLA2)    ,                                             40.55
     &  COMPDA(1,JURSEL)    ,LSWMAT(1,1,JABIN)   ,REFLSO                  40.41 40.03
     &                                                                )   30.21
!TIMG        CALL SWTSTO(117)                                                  40.23
!
!       *** compute transport of action and fill the matrix ***
!
!TIMG        CALL SWTSTA(118)                                                  40.23
        CALL ACTION (IDCMIN      ,IDCMAX            ,SPCSIG            ,  33.09
     &         AC2               ,CAX               ,CAY               ,
     &         CAS               ,CAD               ,SWMATR(1,1,JMATL) ,
     &         SWMATR(1,1,JMATD) ,SWMATR(1,1,JMATU) ,SWMATR(1,1,JMATR) ,
     &         SWMATR(1,1,JMAT5) ,                                        40.22
     &         SWMATR(1,1,JMAT6) ,ISCMIN            ,ISCMAX            ,
     &         IDDLOW            ,IDDTOP            ,ISSTOP            ,
     &                            LSWMAT(1,1,JABLK) ,LSWMAT(1,1,JABIN) ,  30.90
     &         SWMATR(1,1,JLEK1) ,AC1               ,                     40.00
     &         DYNDEP            ,RDX               ,RDY               ,  30.51
     &         SWPDIR            ,IX                ,IY                ,
     &         KSX               ,KSY               ,
     &         XCGRID            ,YCGRID            ,                     40.41
     &         ITER              ,KGRPNT            ,OBREDF            ,  40.41
     &         CAX1              ,CAY1              ,SPCDIR            ,  33.08
     &         CGO               ,SWMATR(1,1,JTRA0) ,SWMATR(1,1,JTRA1)    40.85 33.08
     &                                                                 )
!TIMG        CALL SWTSTO(118)                                                  40.23
!
!       matrix is computed now; updating action densities starts
!       provided ACUPDA is true
!
        IF (.NOT.ACUPDA) THEN                                             40.07
          IF (TESTFL .AND. ITEST.GE.30) WRITE (PRINTF, *) ' No update'    40.07
          GOTO 700                                                        40.07
        ENDIF                                                             40.07
!
!       preparatory steps before solution of linear system
!
!TIMG        CALL SWTSTA(119)                                                  40.23
        CALL SOLPRE(AC2                ,SWMATR(1,1,JAOLD)  ,              40.00
     &              SWMATR(1,1,JMATR)  ,SWMATR(1,1,JMATL)  ,
     &              SWMATR(1,1,JMATD)  ,SWMATR(1,1,JMATU)  ,
     &              SWMATR(1,1,JMAT5)  ,SWMATR(1,1,JMAT6)  ,
     &              IDCMIN             ,IDCMAX             ,
     &              LSWMAT(1,1,JABIN)  ,
     &              IDTOT              ,ISTOT              ,
     &              IDDLOW             ,IDDTOP             ,
     &              ISSTOP             ,
     &              SPCSIG                                 )              40.41 40.23
!TIMG        CALL SWTSTO(119)                                                  40.23
!
        IF ( DYNDEP .OR. ICUR .EQ. 1                                      40.00
     &                                       ) THEN                       30.00
!
!         *** Implicit or explicit scheme in frequency space and ***
!         *** implicit scheme in directional space               ***
!
          IF ( INT(PNUMS(8)) .EQ. 1 ) THEN
!
!           *** Implicit scheme in frequency space. Solve penta- ***
!           *** diagonal system with the SIP solver              ***
!
!TIMG            CALL SWTSTA(120)                                              40.23
            CALL SWSIP ( AC2, SWMATR(1,1,JMATD), SWMATR(1,1,JMATR),
     &                   SWMATR(1,1,JMATL), SWMATR(1,1,JMATU),
     &                   SWMATR(1,1,JMAT5), SWMATR(1,1,JMAT6),
     &                   SWMATR(1,1,JAOLD),
     &                   PNUMS(12), NINT(PNUMS(14)), NINT(PNUMS(13)),
     &                   INOCNV, IDDLOW, IDDTOP, ISSTOP, IDCMIN,          40.41
     &                   IDCMAX )
!TIMG            CALL SWTSTO(120)                                              40.23
!
          ELSE IF (INT(PNUMS(8)).EQ.2 .OR. INT(PNUMS(8)).EQ.3) THEN       40.00
!
!           *** Explicit scheme in frequency space. Energy near the ***
!           *** blocking point is removed from the spectrum based   ***
!           *** on CFL criterion                                    ***
!
!TIMG            CALL SWTSTA(120)                                              40.23
            CALL SOLMT1  (IDCMIN             ,IDCMAX             ,        40.00
     &                    AC2                ,SWMATR(1,1,JMATR)  ,
     &                    SWMATR(1,1,JMATD)  ,SWMATR(1,1,JMATU)  ,
     &                    SWMATR(1,1,JMATL)  ,                            40.23
     &                    ISSTOP             ,                            40.41 40.23
     &                    LSWMAT(1,1,JABLK)  ,IDDLOW             ,        30.90
     &                    IDDTOP                                 )        5/mar
!TIMG            CALL SWTSTO(120)                                              40.23
!
          END IF
!
        ELSE                                                              40.00
!
!         *** No current. Only implicit scheme in directional space  ***
!         *** Solve the tri-diagonal matrix with Thomas algorithm    ***
!
!TIMG          CALL SWTSTA(120)                                                40.23
          CALL SOLMAT (IDCMIN            ,IDCMAX             ,            40.00
     &                AC2                ,SWMATR(1,1,JMATR)  ,
     &                SWMATR(1,1,JMATD)  ,SWMATR(1,1,JMATU)  ,
     &                SWMATR(1,1,JMATL)                                   40.23
     &                                                       )
!TIMG          CALL SWTSTO(120)                                                40.23
!
        END IF
!
!       *** test output ***
!
        IF ( TESTFL .AND. ITEST .GE. 90 ) THEN
          WRITE (PRTEST, *) ' solution vector'                            40.00
          WRITE (PRTEST, *) ' IS ID1 ID2     action densities'            40.00
          DO IS = 1, MSC
            ID_MIN = IDCMIN(IS)
            ID_MAX = IDCMAX(IS)
            WRITE(PRINTF,6621) IS, ID_MIN, ID_MAX,
     &        (AC2(MOD(IDDUM-1+MDC,MDC)+1, IS, KCGRD(1)),                 40.00
     &        IDDUM = ID_MIN, ID_MAX)
 6621       FORMAT(3I4,600(1X,E12.4))                                     40.03
          ENDDO
        END IF
!
!       *** if negative action density occur rescale with a factor ***
!       *** only the sector computed is rescaled !!                ***
!
!TIMG        CALL SWTSTA(121)                                                  40.23
        IF (BRESCL) CALL RESCALE(AC2, ISSTOP, IDCMIN, IDCMAX, NRSCAL)     40.23 40.00
!TIMG        CALL SWTSTO(121)                                                  40.23
!
!       calculate propagation, generation, dissipation, redistribution
!       leak and radiation stress in present grid point
!
!TIMG        CALL SWTSTA(124)                                                  40.23
        IF ( LADDS )                                                      40.85
     &     CALL ADDDIS (COMPDA(1,JDISS)    ,COMPDA(1,JLEAK)    ,
     &                  AC2                ,LSWMAT(1,1,JABIN)  ,          40.02
     &                  SWMATR(1,1,JDIS0)  ,SWMATR(1,1,JDIS1)  ,
     &                  SWMATR(1,1,JGEN0)  ,SWMATR(1,1,JGEN1)  ,          40.85
     &                  SWMATR(1,1,JRED0)  ,SWMATR(1,1,JRED1)  ,          40.85
     &                  SWMATR(1,1,JTRA0)  ,SWMATR(1,1,JTRA1)  ,          40.85
     &                  SWMATR(1,1,JMATL)  ,SWMATR(1,1,JMATU)  ,          40.85
     &                  SWMATR(1,1,JMAT5)  ,SWMATR(1,1,JMAT6)  ,          40.85
     &                  COMPDA(1,JDSXB)    ,                              40.67 40.61
     &                  COMPDA(1,JDSXS)    ,                              40.67 40.61
     &                  COMPDA(1,JDSXW)    ,                              40.67 40.61
     &                  COMPDA(1,JDSXV)    ,                              40.67 40.61
     &                  COMPDA(1,JGSXW)    ,COMPDA(1,JGENR)    ,          40.85
     &                  COMPDA(1,JRSXQ)    ,COMPDA(1,JRSXT)    ,          40.85
     &                  COMPDA(1,JREDS)    ,                              40.85
     &                  COMPDA(1,JTSXG)    ,COMPDA(1,JTSXT)    ,          40.85
     &                  COMPDA(1,JTSXS)    ,COMPDA(1,JTRAN)    ,          40.85
     &                  SWMATR(1,1,JLEK1)  ,COMPDA(1,JRADS)    ,          40.85
     &                  SPCSIG                                            30.72
     &                                                         )
!TIMG        CALL SWTSTO(124)                                                  40.23
!
!       limit the change of the spectrum
!
!TIMG        CALL SWTSTA(122)                                                  40.23
        IF (PNUMS(20).LT.100.) THEN
           IF (IWIND.NE.4 .OR. NSTATC.NE.1) THEN                          40.96 40.61
!             default limiter
              CALL PHILIM (AC2, SWMATR(1,1,JAOLD),                        40.31 40.16
     &                     CGO, KWAVE,                                    40.00
     &                     SPCSIG, LSWMAT(1,1,JABIN),
     &                     ISLMIN, NFLIM,                                 40.23
     &                     QBLOC)                                         30.82
           ELSE
!             Hersbach and Janssen (1999) limiter                         40.61
              CALL HJLIM (AC2, SWMATR(1,1,JAOLD),                         40.61
     &                    CGO, KWAVE,                                     40.61
     &                    SPCSIG, LSWMAT(1,1,JABIN),                      40.61
     &                    ISLMIN, NFLIM,                                  40.61
     &                    QBLOC, COMPDA(1,JUSTAR))                        40.61
           END IF
        END IF
!TIMG        CALL SWTSTO(122)                                                  40.23
!
!       *** reduce the computed energy density if the value is  ***
!       *** larger then the limit value as computed in SWIND    ***
!       *** in case of first or second generation mode          ***
!
!TIMG        CALL SWTSTA(123)                                                  40.23
        IF ( IWIND .EQ. 1 .OR. IWIND .EQ. 2 )
     &    CALL WINDP3 (ISSTOP, ALIMW, AC2, GROWW, IDCMIN, IDCMAX )        40.41 30.21
!TIMG        CALL SWTSTO(123)                                                  40.23
!
!       *** test output ***
!
        IF ( TESTFL .AND. ITEST .GE. 70 ) THEN
          WRITE (PRINTF, *) ' action densities after adaptations'         40.00
          WRITE (PRTEST, *) ' IS ID1 ID2     action densities'            40.00
          DO IS = 1, MSC
            ID_MIN = IDCMIN(IS)
            ID_MAX = IDCMAX(IS)
            WRITE(PRINTF,6621) IS, ID_MIN, ID_MAX,
     &        (AC2(MOD(IDDUM-1+MDC,MDC)+1, IS, KCGRD(1)),                 40.00
     &        IDDUM = ID_MIN, ID_MAX)
          ENDDO
        END IF
!
 700    CONTINUE                                                          40.07
!
      END IF
!
!     End of the subroutine SWOMPU
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWPRSET (SPCSIG)
!
!****************************************************************
!
      USE OCPCOMM4
      USE SWCOMM2
      USE SWCOMM3
      USE SWCOMM4
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
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     40.80, Oct. 07: New subroutine
!
!  2. Purpose
!
!     Print all the settings used in SWAN run
!
!  4. Argument variables
!
!     SPCSIG      Relative frequencies in computational domain in sigma-space
!
      REAL SPCSIG(MSC)
!
!  6. Local variables
!
!     ICMX  :     stencil size
!     IENT  :     number of entries in this subroutine
!     IS    :     loop counter
!
      INTEGER ICMX, IS, IENT
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPRSET')

      WRITE(PRINTF,333) 'SWAN'                                            40.00
 333  FORMAT(/,
     &'----------------------------------------------------------------'
     &,/,
     &'                  COMPUTATIONAL PART OF ', A
     &,/,
     &'----------------------------------------------------------------'
     &,/)
!
      IF ( ONED ) THEN                                                    32.02
         WRITE(PRINTF,*) 'One-dimensional mode of SWAN is activated'      32.02
      ENDIF                                                               32.02
!
      IF (OPTG.NE.5) THEN                                                 40.80
         WRITE(PRINTF,7001) MXC,MYC
 7001    FORMAT(' Gridresolution       : MXC    ',I12  ,' MYC   ',I12)
         WRITE(PRINTF,7002) MCGRD                                         40.41
 7002    FORMAT('                      : MCGRD  ',I12)                    40.41
         WRITE(PRINTF,7101) MSC,MDC
 7101    FORMAT('                      : MSC    ',I12  ,' MDC   ',I12)
      ELSE
         WRITE(PRINTF,7102) MSC,MDC
 7102    FORMAT(' Gridresolution       : MSC    ',I12  ,' MDC   ',I12)
      ENDIF
      WRITE(PRINTF,7201) MTC                                              40.41
 7201 FORMAT('                      : MTC    ',I12)                       40.41
      WRITE(PRINTF,7301) NSTATC, ITERMX
 7301 FORMAT('                      : NSTATC ',I12  ,' ITERMX',I12)
      WRITE(PRINTF,7013) ITFRE,IREFR
 7013 FORMAT(' Propagation flags    : ITFRE  ',I12  ,' IREFR ',I12)
      WRITE(PRINTF,7014) IBOT,ISURF
 7014 FORMAT(' Source term flags    : IBOT   ',I12  ,' ISURF ',I12)
      WRITE(PRINTF,7114) IWCAP,IWIND
 7114 FORMAT('                      : IWCAP  ',I12  ,' IWIND ',I12)
      WRITE(PRINTF,7015) ITRIAD,IQUAD
 7015 FORMAT('                      : ITRIAD ',I12  ,' IQUAD ',I12)
      WRITE(PRINTF,7017) IVEG                                             40.55
 7017 FORMAT('                      : IVEG   ',I12)                       40.55
      IF (OPTG.NE.5) WRITE(PRINTF,7004) DX,DY                             40.80
 7004 FORMAT(' Spatial step         : DX     ',E12.4,' DY    ',E12.4)
      WRITE(PRINTF,7104) EXP(ALOG(SHIG/SLOW)/REAL(MSC-1))-1.,             40.41
     &                   DDIR/DEGRAD                                      40.41
 7104 FORMAT(' Spectral bin         : df/f   ',E12.4,' DDIR  ',E12.4)     40.41
      WRITE(PRINTF,7003) GRAV, RHO                                        40.41
 7003 FORMAT(' Physical constants   : GRAV   ',E12.4,' RHO   ',E12.4)     40.41
      WRITE(PRINTF,7027) U10 , WDIC/DEGRAD
 7027 FORMAT(' Wind input           : WSPEED ',E12.4,' DIR   ',E12.4)
      WRITE(PRINTF,7123) PWTAIL(1),PWTAIL(2)
 7123 FORMAT(' Tail parameters      : E(f)   ',E12.4,' E(k)  ',E12.4)
      WRITE(PRINTF,7133) PWTAIL(3),PWTAIL(4)
 7133 FORMAT('                      : A(f)   ',E12.4,' A(k)  ',E12.4)
      WRITE(PRINTF,8013) PNUMS(1), PNUMS(4)                               40.41
 8013 FORMAT(' Accuracy parameters  : DREL   ',E12.4,' NPNTS ',E12.4)     40.41
      IF (PNUMS(21).EQ.0.) THEN                                           40.41
         WRITE(PRINTF,8213) PNUMS(15),PNUMS(16)                           40.41
      ELSE IF (PNUMS(21).EQ.1.) THEN                                      40.41
         WRITE(PRINTF,8214) PNUMS(2),PNUMS(15)                            40.41
      END IF                                                              40.41
 8213 FORMAT('                      : DHOVAL ',E12.4,' DTOVAL',E12.4)     40.41
 8214 FORMAT('                      : DHABS  ',E12.4,' CURVAT',E12.4)     40.41
      WRITE(PRINTF,3613) PNUMS(20)
 3613 FORMAT('                      : GRWMX  ',E12.4)
      WRITE(PRINTF,8508) WLEV, DEPMIN                                     40.41
 8508 FORMAT(' Drying/flooding      : LEVEL  ',E12.4,' DEPMIN',E12.4)     40.41
      IF (BNAUT) THEN                                                     40.41
         WRITE (PRINTF,8510) 'nautical '                                  40.41
      ELSE                                                                40.41
         WRITE (PRINTF,8510) 'Cartesian'                                  40.41
      ENDIF                                                               40.41
 8510 FORMAT(' The ',A9,                                                  40.41
     &       ' convention for wind and wave directions is used')          40.41
      IF (OPTG.EQ.5) THEN                                                 40.80
        WRITE(PRINTF,8511) 'BSBT  '                                       40.80
      ELSEIF (PROPSC.EQ.3) THEN                                           40.80 40.41
        WRITE(PRINTF,8511) 'S&L   '                                       40.41
        ICMX = 10                                                         40.41
      ELSEIF (PROPSC.EQ.2) THEN                                           40.41
        WRITE(PRINTF,8511) 'SORDUP'                                       40.41
        ICMX = 5                                                          40.41
      ELSE                                                                40.41
        WRITE(PRINTF,8511) 'BSBT  '                                       40.41
        ICMX = 3                                                          40.41
      ENDIF                                                               40.41
 8511 FORMAT(' Scheme for geographic propagation is ',A6)                 40.41
      IF (OPTG.NE.5) WRITE(PRINTF,8512) PROPSC, ICMX                      40.80 40.41
 8512 FORMAT(' Scheme geogr. space  : PROPSC ',I12  ,' ICMAX ',I12)       40.41
      WRITE(PRINTF,8513) PNUMS(7), PNUMS(6)                               40.41
 8513 FORMAT(' Scheme spectral space: CSS    ',E12.4,' CDD   ',E12.4)     40.41
!
      IF ( (DYNDEP .OR. ICUR.EQ.1) .AND. INT(PNUMS(8)).EQ.1 ) THEN        40.41
         WRITE(PRINTF,*) 'Solver is SIP'                                  40.23
         WRITE(PRINTF,2213) PNUMS(12), INT(PNUMS(13))
 2213    FORMAT('                      : EPS2   ',E12.4,' OUTPUT',I12)
         WRITE(PRINTF,8223) INT(PNUMS(14))
 8223    FORMAT('                      : NITER  ',I12)
      ENDIF
!
      IF (ICUR.GT.0) THEN
         WRITE (PRINTF,7115) 'on'
      ELSE
         WRITE (PRINTF,7115) 'off'
      ENDIF
 7115 FORMAT(' Current is ', A3)
!
      IF (IQUAD.GT.0) THEN
         WRITE(PRINTF,8413) IQUAD
 8413    FORMAT(' Quadruplets          : IQUAD  ',I12)
         IF (IQUAD.LE.3 .OR. IQUAD.EQ.8) THEN                             40.41
            WRITE(PRINTF,8414) PQUAD(1), PQUAD(2)                         40.41
 8414       FORMAT('                      : LAMBDA ',E12.4,               40.41
     &                                    ' CNL4  ',E12.4)                40.41
            WRITE(PRINTF,8415) PQUAD(3), PQUAD(4)                         40.41
 8415       FORMAT('                      : CSH1   ',E12.4,               40.41
     &                                    ' CSH2  ',E12.4)                40.41
            WRITE(PRINTF,8416) PQUAD(5)                                   40.41
 8416       FORMAT('                      : CSH3   ',E12.4)               40.41
         ENDIF                                                            40.41
         WRITE(PRINTF,8417) PTRIAD(3)                                     40.41
 8417    FORMAT(' Maximum Ursell nr for Snl4 :  ',E12.4)                  40.41
      ELSE
         WRITE (PRINTF, *) 'Quadruplets is off'
      ENDIF                                                               40.41
      IF (ITRIAD.GT.0) THEN                                               40.41
         WRITE(PRINTF,9413) ITRIAD, PTRIAD(1)
 9413    FORMAT(' Triads               : ITRIAD ',I12  ,
     &                                 ' TRFAC ',E12.4)
         WRITE(PRINTF,9414) PTRIAD(2), PTRIAD(4)                          40.41
 9414    FORMAT('                      : CUTFR  ',E12.4,                  40.41
     &                                 ' URCRI ',E12.4)                   40.41
         WRITE(PRINTF,9415) PTRIAD(5)                                     40.41
 9415    FORMAT(' Minimum Ursell nr for Snl3 :  ',E12.4)                  40.41
      ELSE
         WRITE (PRINTF, *) 'Triads is off'
      ENDIF                                                               40.41
      IF (IBOT.EQ.2) THEN
         WRITE(PRINTF,7005) PBOT(2), PBOT(1)
 7005    FORMAT(' Collins (`72)        : CFW    ',E12.4,' CFC   ',E12.4)
      ELSEIF (IBOT.EQ.3) THEN
         WRITE(PRINTF,7335) PBOT(4), PBOT(5)
 7335    FORMAT(' Madsen et al. (`84)  : MF     ',E12.4,' KN    ',E12.4)
      ELSEIF (IBOT.EQ.1) THEN
         WRITE(PRINTF,7325) PBOT(3)
 7325    FORMAT(' JONSWAP (`73)        : GAMMA  ',E12.4)
      ELSEIF (IBOT.EQ.4) THEN                                             41.04
         WRITE(PRINTF,7326) PBOT(6), PBOT(7)
 7326    FORMAT(' JONSWAP (`73)        : GAMMA1 ',E12.4,' GAMMA2',E12.4)
         WRITE(PRINTF,7327) PBOT(8), PBOT(9)
 7327    FORMAT('                      : DSPR1  ',E12.4,' DSPR2 ',E12.4)
      ELSE
         WRITE (PRINTF, *) 'Bottom friction is off'
      ENDIF
!
      IF (IVEG.EQ.1) THEN                                                 40.55
        WRITE(PRINTF,7008)                                                40.55
 7008   FORMAT(' Vegetation due to Dalrymple (1984)')                     40.55
      ELSE                                                                40.55
        WRITE (PRINTF, *) 'Vegetation is off'                             40.55
      ENDIF                                                               40.55
!
      IF (IWCAP.EQ.1) THEN
         WRITE(PRINTF,6005) PWCAP(1),PWCAP(2)
 6005    FORMAT(' W-cap Komen (`84)    : EMPCOF ',E12.4,
     &          ' APM   ',E12.4)
      ELSEIF (IWCAP.EQ.2) THEN
         WRITE(PRINTF,6335) PWCAP(3),PWCAP(4)
 6335    FORMAT(' W-cap Janssen (`90)  : CFJANS ',E12.4,
     &          ' DELTA ',E12.4)
      ELSEIF (IWCAP.EQ.3) THEN
         WRITE(PRINTF,6135) PWCAP(5)
 6135    FORMAT(' W-cap Longuet-Higgins: CFLHIG ',E12.4)
      ELSEIF (IWCAP.EQ.4) THEN
         WRITE(PRINTF,6136) PWCAP(6), PWCAP(7)
 6136    FORMAT(' W-cap Battjes/Janssen: BJSTP  ',E12.4,
     &          ' BJALF ',E12.4)
      ELSEIF (IWCAP.EQ.5) THEN
         WRITE(PRINTF,6136) PWCAP(6), PWCAP(7)
         WRITE(PRINTF,6137) PWCAP(8)
 6137    FORMAT('                      : KCONV  ',E12.4)
      ELSEIF (IWCAP.EQ.7) THEN                                            40.53
         WRITE(PRINTF,6139) PWCAP(1), PWCAP(12)                           40.53
         WRITE(PRINTF,6140) PWCAP(9), PWCAP(11)                           40.53
 6139    FORMAT(' W-cap Alves-Banner   : CDS2   ',                        40.53
     &          E12.4,' BR    ',E12.4)                                    40.53
 6140    FORMAT('                      : POWST  ',                        40.53
     &          E12.4,' POWK  ',E12.4)                                    40.53
      ELSE
         WRITE (PRINTF, *) 'Whitecapping is off'
      ENDIF
!
      IF (ISURF.EQ.1) THEN
         WRITE(PRINTF,7012) PSURF(1),PSURF(2)
 7012    FORMAT(' Battjes&Janssen (`78): ALPHA  ',E12.4,
     &          ' GAMMA ',E12.4)
      ELSEIF (ISURF.EQ.2) THEN                                            970219
         WRITE(PRINTF,7212) PSURF(1), PSURF(4)
 7212    FORMAT(' Nelson (`94)         : ALPHA  ',E12.4,
     &          ' GAMmin',E12.4)
         WRITE(PRINTF,7213) PSURF(5)
 7213    FORMAT('                        GAMmax ',E12.4)
      ELSEIF (ISURF.EQ.3) THEN                                            41.03
         WRITE(PRINTF,7214) PSURF(1), PSURF(4)
 7214    FORMAT(' Ruessink et al (2003): ALPHA  ',E12.4,
     &          ' A     ',E12.4)
         WRITE(PRINTF,7215) PSURF(5)
 7215    FORMAT('                        B      ',E12.4)
      ELSEIF (ISURF.EQ.4) THEN                                            41.03
         WRITE(PRINTF,7216) PSURF(1), PSURF(2)
 7216    FORMAT(' Thornton&Guza (`83)  : ALPHA  ',E12.4,
     &          ' GAMMA ',E12.4)
         WRITE(PRINTF,7217) PSURF(4)
 7217    FORMAT('                        N      ',E12.4)
      ELSE
         WRITE (PRINTF, *) 'Surf breaking is off'
      ENDIF
!
      IF (LSETUP.GT.0) THEN
         WRITE(PRINTF,7109) PSETUP(2)
 7109    FORMAT(' Set-up               : SUPCOR ',E12.4)
      ELSE
         WRITE (PRINTF, *) 'Set-up is off'
      ENDIF
!
      IF (IDIFFR.EQ.1) THEN                                               40.41
        WRITE(PRINTF,7110) PDIFFR(1), NINT(PDIFFR(2))
 7110   FORMAT(' Diffraction          : SMPAR  ',E12.4,
     &         ' SMNUM ',I12)
      ELSE
        WRITE (PRINTF, *) 'Diffraction is off'
      ENDIF
!
      WRITE(PRINTF,7126) PWIND(14), PWIND(15)
 7126 FORMAT(' Janssen (`89,`90)    : ALPHA  ',E12.4,
     &       ' KAPPA ',E12.4)
      WRITE(PRINTF,7136) PWIND(16), PWIND(17)
 7136 FORMAT(' Janssen (`89,`90)    : RHOA   ',E12.4,
     &       ' RHOW  ',E12.4)
      WRITE(PRINTF,*)
      WRITE(PRINTF,1012) PWIND(1), PWIND(2)
 1012 FORMAT(' 1st and 2nd gen. wind: CF10   ',E12.4,
     &       ' CF20  ',E12.4)
      WRITE(PRINTF,1013) PWIND(3), PWIND(4)
 1013 FORMAT('                      : CF30   ',E12.4,
     &       ' CF40  ',E12.4)
      WRITE(PRINTF,1014) PWIND(5), PWIND(6)
 1014 FORMAT('                      : CF50   ',E12.4,
     &       ' CF60  ',E12.4)
      WRITE(PRINTF,1015) PWIND(7), PWIND(8)
 1015 FORMAT('                      : CF70   ',E12.4,
     &       ' CF80  ',E12.4)
      WRITE(PRINTF,1016) PWIND(9), PWIND(10)
 1016 FORMAT('                      : RHOAW  ',E12.4,
     &       ' EDMLPM',E12.4)
      WRITE(PRINTF,1017) PWIND(11), PWIND(12)
 1017 FORMAT('                      : CDRAG  ',E12.4,
     &       ' UMIN  ',E12.4)
      WRITE(PRINTF,1018) PWIND(13)
 1018 FORMAT('                      : LIM_PM ',E12.4)
!
      WRITE(PRINTF,*)
      IF ( ITEST .GT. 2 )  THEN
         DO IS = 1, MSC
            WRITE(PRINTF,*)' IS and SPCSIG(IS)    :',IS,SPCSIG(IS)        30.72
         ENDDO
         WRITE(PRINTF,*)
      ENDIF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SACCUR (DEP2       ,KGRPNT     ,
     &                   XYTST      ,                                     40.41
     &                   AC2        ,SPCSIG     ,ACCUR      ,             30.72
     &                   HSACC1     ,HSACC2     ,SACC1      ,
     &                   SACC2      ,DELHS      ,DELTM      ,             40.31
     &                   I1MYC      ,I2MYC                  )             40.31
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
      USE M_PARALL                                                        40.31
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer(s): R.C. Ris                                   |
!     |                Modified by R. Padilla and N. Booij        |
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
!     40.03: Nico Booij
!     40.22: John Cazes and Tim Campbell
!     40.30: Marcel Zijlema
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!
!  1. Update
!
!     30.72, Nov. 97: Declaration of DDIR, PI and PI2 removed because
!                     they are common and already declared in the
!                     INCLUDE file
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Aug. 99: Introduced a new overall measure for checking accuracy
!     30.82, Aug. 99: Changed all variables INDEX to INDX, since INDEX is reserved
!     40.03, Feb. 00: test level of message changed
!     40.22, Sep. 01: Added initialization of SACC1 and HSACC1 elements   40.22
!                     that are not wet points.                            40.22
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Jul. 03: some improvements and corrections w.r.t. OpenMP
!     40.41, Aug. 04: add some test output for checking accuracy
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     To check the accuracy of the final computation. If a certain
!     accuracy has been reached then terminate the iteration process
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     AC2         action density
!     ACCUR       indicates percentage of grid points in
!                 which accuracy is reached
!     DELHS       difference in Hs between last 2 iterations
!     DELTM       difference in Tm01 between last 2 iterations
!     DEP2        depth
!     HSACC1      significant wave height at iter-1
!     HSACC2      significant wave height at iter
!     I1MYC       lower index for thread loop over y-grid row
!     I2MYC       upper index for thread loop over y-grid row
!     KGRPNT      indirect addressing
!     SACC1       mean wave frequency at iter-1
!     SACC2       mean wave frequency at iter
!     SPCSIG      relative frequencies in computational domain
!                 in sigma-space                                          30.72
!     XYTST       coordinates of test points
!
      INTEGER I1MYC, I2MYC                                                40.31
      REAL    ACCUR
      INTEGER KGRPNT(MXC,MYC)
      INTEGER XYTST(2*NPTST)                                              40.41
      REAL    AC2(MDC,MSC,MCGRD)   ,
     &        DEP2(MCGRD)          ,
     &        HSACC1(MCGRD)        ,
     &        HSACC2(MCGRD)        ,
     &        SACC1(MCGRD)         ,
     &        SACC2(MCGRD)         ,
     &        DELHS(MCGRD)         ,
     &        DELTM(MCGRD)
      REAL    SPCSIG(MSC)                                                 30.72
!
!  6. Local variables
!
      INTEGER  IS    ,ID    ,WETGRD,IACCUR,IX,IY,IX1,IX2,IY1,IY2
      INTEGER  WETGRDt, IACCURt                                           40.31
!
!     INDX  : counter                                                     30.82
!     NINDX : number of gridpoints to average over                        30.82
!
      INTEGER INDX, NINDX                                                 30.82
      INTEGER NINDXt                                                      40.31
!
      REAL    SME_T ,SME_B ,                                              30.72
     &        TMREL ,HSREL ,TMABS ,HSABS                                  30.72
      REAL    ARR(2)                                                      40.30
!
!     HSMN2 : mean Hs over space at current iteration level               30.82
!     HSOVAL: Overall accuracy measure for Hs                             30.82
!     SMN2  : mean Tm over space at current iteration level               30.82
!     TMOVAL: Overall accuracy measure for Tm                             30.82
!
      REAL    HSMN2, HSOVAL, SMN2, TMOVAL                                 30.82
      REAL    HSMN2t, SMN2t                                               40.31
!
!     LHEAD : logical indicating to write header                          40.41
!     TSTFL : indicates whether grid point is a test point                40.41
!
      LOGICAL LHEAD, TSTFL                                                40.41
!
!  7. Common blocks used
!
!
!     Place local summed variables in common block so they will           40.31
!     be scoped as shared                                                 40.31
      COMMON/SACCUR_MT_COM/HSMN2,SMN2,NINDX,WETGRD,IACCUR                 40.31
!
!  8. Subroutines used
!
!     EQREAL           Boolean function which compares two REAL values
!     STRACE           Tracing routine for debugging
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     SWREDUCE         Performs a global reduction
!
      LOGICAL EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP (in SWANCOM1)
!
! 12. Structure
!
!   ---------------------------------------------------------------
!   If not the first iteration, the do
!     Set old values in dummy array
!   ---------------------------------------------------------------
!   Do for every x and y
!     Compute the mean action density frequency SACC1 and the
!     and the significant waveheight HSACC1
!   ---------------------------------------------------------------
!   If relative error for mean frequency or significant wave height
!      > certain given value then increase variable with one and
!      compute the relative number of gridpoints in where the accuracy
!      has not been reached
!   ---------------------------------------------------------------
!   End of the subroutine SACCUR
!   ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SACCUR')
!
!$OMP MASTER                                                              40.31
!     Master thread initialize the shared variables                       40.31
      HSMN2  = 0.                                                         40.31
       SMN2  = 0.                                                         40.31
      NINDX  = 0                                                          40.31
      WETGRD = 0                                                          40.31
      IACCUR = 0                                                          40.31
!$OMP END MASTER                                                          40.31
!$OMP BARRIER                                                             40.31
!
      IF ( LMXF ) THEN                                                    40.41 40.30
         IX1 = 1
      ELSE
         IX1 = 1+IHALOX
      END IF
      IF ( LMXL ) THEN
         IX2 = MXC
      ELSE
         IX2 = MXC-IHALOX
      END IF
      IF ( LMYF ) THEN
         IY1 = I1MYC
      ELSE
         IY1 = 1+IHALOY
      END IF
      IF ( LMYL ) THEN
         IY2 = I2MYC
      ELSE
         IY2 = MYC-IHALOY
      END IF
!
!     *** If the computation is non steady : check the gridpoints ***
!     *** at the different "timesteps" if they are still the same ***
!     *** then WETGRD is the same. If not: change this subroutine ***
!     ***                                                         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***   ++++++++     ++++++     +++++++++      ++++++         ***
!     ***   +                 +     ++++                +         ***
!     ***   +                 +     +++                 +         ***
!     ***   +                 +     ++                  +         ***
!     ***   +                 +     +                   +         ***
!     ***   +       t=0       +     +      t=t+1        +         ***
!     ***   +                 +     +                   +         ***
!     ***   +++++++++++++++++++     +++++++++++++++++++++         ***
!     ***                                                         ***
!
      WETGRDt = 0                                                         40.31
      DO 101 IX = IX1, IX2
      DO 100 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
        IF (DEP2(INDX) .GT. DEPMIN) THEN                                  25/MAR
          HSACC1(INDX) = MAX( 1.E-20 , HSACC2(INDX) )                     30.21
          SACC1(INDX)  = MAX( 1.E-20 , SACC2(INDX)  )                     30.21
          WETGRDt = WETGRDt + 1                                           40.31
!       Added to initialize HSACC1 and SACC1 values.                      40.22
        ELSE                                                              40.22
          HSACC1(INDX) = 0.                                               40.22
          SACC1(INDX)  = 0.                                               40.22
        END IF                                                            40.22
 100  CONTINUE
 101  CONTINUE
!
!     *** first criterion to terminate the iteration process ***
!     ***                                                    ***
!     *** RELATIVE error :                                   ***
!     ***               Hs2  - Hs1      Tm2 - Tm1            ***
!     ***      DREL  =  ----------  and ----------           ***
!     ***                   Hs1            Tm1               ***
!     ***                                                    ***
!     *** ABSOLUTE error :                                   ***
!     ***                                                    ***
!     ***      DHABS  =  Hs2  - Hs1 < PNUMS(2)               ***
!     ***      DTABS  =  Tm2  - Tm1 < PNUMS(3)               ***
!     ***                                                    ***
!
      DO 201 IX = IX1, IX2
      DO 200 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
!
!       *** Compute the mean ENERGY DENSITY frequency and    ***
!       *** significant waveheight over the full spectrum    ***
!       *** per gridpoint                                    ***
!
        IF (DEP2(INDX) .GT. DEPMIN) THEN                                  25/MAR
          SME_T  = 0.
          SME_B  = 0.
          DO 180 IS = 1, MSC
            DO 170 ID = 1, MDC
              ACS2  = SPCSIG(IS)**2 * AC2(ID,IS,INDX)                     30.72
              ACS3  = SPCSIG(IS) * ACS2                                   30.72
              SME_B = SME_B + ACS2
              SME_T = SME_T + ACS3
 170        CONTINUE
 180      CONTINUE
          SME_B = SME_B * FRINTF * DDIR                                   30.50
          SME_T = SME_T * FRINTF * DDIR                                   30.50
!
!         *** mean frequency and significant wave height per gridpoint ***
!
          IF ( SME_B .LE. 0. ) THEN
            SME_B = 1.E-20
            SACC2(INDX) = 1.E-20                                          30.21
            HSACC2(INDX) = 1.E-20                                         30.21
          ELSE
            SACC2(INDX) = MAX ( 1.E-20 , (SME_T / SME_B) )                30.21
            HSACC2(INDX) = MAX ( 1.E-20 , (4. * SQRT(SME_B)) )            30.21
          END IF
        END IF
 200  CONTINUE
 201  CONTINUE
!
!     *** the mean significant waveheight and the mean  ***
!     *** relative frequency over the gridpoints which  ***
!     *** depth is larger than 0 m.                     ***
!     *** The amount of gridpoints is denoted with the  ***
!     *** variable : NINDX                              ***
!     *** These values are used to compute the SRELF    ***
!     *** and the HSRELF instead of SACC2 and HSACC2    ***
!
!     Note that initialization of ACCUR is not needed since it is         40.31
!     not being summed upon                                               40.31
      IACCURt = 0                                                         40.31
!
!     Calculate the mean Hs and Tm over all wet gridpoints. These means
!     are then used as an overall accuracy measure.
!
      HSMN2t= 0.                                                          40.31 30.82
       SMN2t= 0.                                                          40.31 30.82
      NINDXt= 0                                                           40.31 30.82
!
      DO 301 IX = IX1, IX2
      DO 300 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)
        IF (DEP2(INDX).GT.DEPMIN) THEN                                    30.82
          HSMN2t = HSMN2t + HSACC2(INDX)                                  40.31 30.82
          SMN2t  =  SMN2t +  SACC2(INDX)                                  40.31 30.82
          NINDXt = NINDXt + 1                                             40.31 30.82
        END IF                                                            30.82
 300  CONTINUE                                                            30.82
 301  CONTINUE
!
!     Global sum of NINDX                                                 40.31
!$OMP ATOMIC                                                              40.31
      NINDX = NINDX + NINDXt                                              40.31
!$OMP ATOMIC                                                              40.31
      HSMN2 = HSMN2 + HSMN2t                                              40.31
!$OMP ATOMIC                                                              40.31
      SMN2 = SMN2 + SMN2t                                                 40.31
!
!$OMP BARRIER                                                             40.31
!$OMP MASTER                                                              40.31
      ARR(1) = HSMN2                                                      40.30
      ARR(2) = SMN2                                                       40.30
      CALL SWREDUCE( ARR, 2, SWREAL, SWSUM )                              40.30
      CALL SWREDUCE( NINDX, 1, SWINT, SWSUM )                             40.41
      IF (STPNOW()) RETURN                                                40.30
      HSMN2 = ARR(1) / REAL(NINDX)                                        40.30 30.82
       SMN2 = ARR(2) / REAL(NINDX)                                        40.30 30.82
!$OMP END MASTER                                                          40.31
!$OMP BARRIER                                                             40.31
!
!     Calculate a set of accuracy parameters based on relative, absolute
!     and overall accuracy measures for Hs and Tm
!
      LHEAD=.TRUE.                                                        40.41
      DO 401 IX = IX1, IX2
      DO 400 IY = IY1, IY2
        INDX = KGRPNT(IX,IY)

!       --- determine whether the point is a test point                   40.41

        TSTFL = .FALSE.                                                   40.41
        IF (NPTST.GT.0) THEN                                              40.41
           DO 20 II = 1, NPTST                                            40.41
              IF (IX.NE.XYTST(2*II-1)) GOTO 20                            40.41
              IF (IY.NE.XYTST(2*II  )) GOTO 20                            40.41
              TSTFL = .TRUE.                                              40.41
  20       CONTINUE                                                       40.41
        END IF                                                            40.41

        IF ( DEP2(INDX) .GT. DEPMIN ) THEN                                25/MAR
          TMREL  = ABS ( SACC2(INDX) - SACC1(INDX) ) /
     &                   SACC1(INDX)                                      30.21
          TMABS  = ABS ( ( PI2/SACC2(INDX)) - (PI2/SACC1(INDX)) )         30.21
          TMOVAL = ABS ( SACC2(INDX) - SACC1(INDX) ) / SMN2               30.82
!
          HSREL  = ABS ( HSACC2(INDX) - HSACC1(INDX) ) /
     &                   HSACC1(INDX)                                     30.21
          HSABS  = ABS ( HSACC2(INDX) - HSACC1(INDX) )                    30.21
          HSOVAL = ABS ( HSACC2(INDX) - HSACC1(INDX) ) / HSMN2            30.82
!
          IF (EQREAL(SACC1(INDX),1.E-20) .OR.                             40.41
     &        EQREAL(SACC2(INDX),1.E-20) ) THEN                           40.41
             DELTM(INDX) = 0.                                             40.41
          ELSE                                                            40.41
             DELTM(INDX) = TMABS
          END IF                                                          40.41
          DELHS(INDX) = HSABS
!
!         *** gridpoint in which mean period and wave height ***
!         *** have reached required accuracy                 ***
!
          IF ( ITEST .GE. 30 .AND. TESTFL) THEN
            WRITE(PRINTF,3002) SACC2(INDX), SACC1(INDX),
     &                         HSACC2(INDX), HSACC1(INDX)                 30.21
 3002       FORMAT(' SACCUR: SA2 SA1 HSA2 HSA1       :',4E12.4)
            WRITE(PRINTF,2002) TMREL, HSREL, TMABS, HSABS
 2002       FORMAT(' SACCUR: TMREL HSREL TMABS HSABS :',4E12.4)
          ENDIF
!
          IF ( (TMREL .LE. PNUMS(1) .OR. TMOVAL .LE. PNUMS(16)) .AND.     30.82
     &         (HSREL .LE. PNUMS(1) .OR. HSOVAL .LE. PNUMS(15)) ) THEN    30.82
            IACCURt = IACCURt + 1                                         40.31
          END IF
!
          IF (TSTFL) THEN                                                 40.41
             IF (LHEAD) WRITE(PRINTF,501)                                 40.41
             WRITE(PRINTF,502) IX+MXF-2, IY+MYF-2, HSREL, HSOVAL,         40.41
     &                         TMREL, TMOVAL                              40.41
 501         FORMAT(25X,'dHrel          ','dHoval         ',              40.41
     &                  'dTm01rel       ','dTm01oval ')                   40.41
 502         FORMAT(1X,SS,'(IX,IY)=(',I5,',',I5,')','  ',                 40.41
     &              1PE13.6E2,'  ',1PE13.6E2,'  ',1PE13.6E2,'  ',         40.41
     &              1PE13.6E2)                                            40.41
             LHEAD=.FALSE.                                                40.41
          END IF                                                          40.41
        ELSE
!         *** otherwise set arrays equal 0 ***
          DELTM(INDX) = 0.0
          DELHS(INDX) = 0.0
        END IF
!
!       Test output at test points
!
        IF ( ITEST .GE. 30 .AND. TESTFL) THEN                             30.82
           WRITE(PRINTF,1003) TMREL, TMABS, TMOVAL                        30.82
 1003      FORMAT(' SACCUR: TMREL, TMABS, TMOVAL  :',3E12.4)              30.82
           WRITE(PRINTF,1004) HSREL, HSABS, HSOVAL                        30.82
 1004      FORMAT(' SACCUR: HSREL, HSABS, HSOVAL  :',3E12.4)              30.82
        END IF                                                            30.82
!
 400  CONTINUE
 401  CONTINUE
!
!     Global sum of IACCUR and WETGRD                                     40.31
!$OMP ATOMIC                                                              40.31
      IACCUR = IACCUR + IACCURt                                           40.31
!$OMP ATOMIC                                                              40.31
      WETGRD = WETGRD + WETGRDt                                           40.31
!
!$OMP BARRIER                                                             40.31
!$OMP MASTER                                                              40.31
!
      CALL SWREDUCE ( IACCUR, 1, SWINT, SWSUM )                           40.30
      IF (STPNOW()) RETURN                                                40.30
      ACCUR  = REAL(IACCUR) * 100. / REAL(NINDX)
!$OMP END MASTER                                                          40.31
!$OMP BARRIER                                                             40.31
!
!     *** test output ***
!
!$OMP MASTER                                                              40.31
      IF ( ITEST .GE. 30 ) THEN                                           40.03
        WRITE(PRINTF,1002) PNUMS(1), PNUMS(2), PNUMS(3)
 1002   FORMAT(' SACCUR: PNUMS(1) DHABS DTABS  :',3E12.4)
        WRITE(PRINTF,1008) NINDX,IACCUR,ACCUR
 1008   FORMAT(' SACCUR: WETGRD IACCUR ACCUR   :',2I8,E12.4)
      END IF
!$OMP END MASTER                                                          40.31
!
!     End of the subroutine SACCUR
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE INSAC (AC2      ,SPCSIG   ,DEP2     ,                    30.72
     &                  HSACC2   ,SACC2    ,KGRPNT   ,                    40.31 40.30
     &                  I1MYC    ,I2MYC              )                    40.31
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
      USE M_PARALL                                                        40.31
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer(s): R.C. Ris                                   |
!     |                Modified by R. Padilla and N. Booij        |
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
!     40.30: Marcel Zijlema
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!
!  1. Update
!
!     30.72, Nov. 97: Declartion of DDIR removed because it is a common
!                     and already declared in the INCLUDE file
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Jul. 03: some improvements and corrections w.r.t. OpenMP
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     To check the accuracy of the final computation. If a certain
!     accuracy has been reached then quit the iteration process
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!        IX          Counter of gridpoints in x-direction
!        IY          Counter of gridpoints in y-direction
!        IS          Counter of relative frequency band
!        ID          Counter of the spectral direction
!        ICMAX       Maximum counter for the points of the molecul
!        ITER        Number of iteration i.e. number of full sweeps
!        MXC         Maximum counter of gridppoints in x-direction
!        MYC         Maximum counter of gridppoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!
!        DDIR        Spectral direction band width
!        DS          Width of the frequency band
!
!        one and more dimensional arrays:
!        ---------------------------------
!
!        AC2       4D    Action density as function of D,S,X,Y at time T
!        DEP2      2D    Depth
!        HSACC2    2D    Dummy array for the significant wave height
!                        (old value)
!        SACC2     2D    Dummy array for the mean frequency (old value)
!
!     5. SUBROUTINES CALLING
!
!        SWOMPU
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
!     9. STRUCTURE
!
!   ---------------------------------------------------------------
!   If the first iteration, the do
!     Set old values in dummy array
!   ---------------------------------------------------------------
!   End of the subroutine INSAC
!   ----------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  IS    ,ID, IX, IY, IX1, IX2, IY1, IY2
      INTEGER  I1MYC, I2MYC                                               40.31
!
      INTEGER  KGRPNT(MXC,MYC)
      REAL     SME_T ,SME_B                                               30.72
!
      REAL     AC2(MDC,MSC,MCGRD)   ,                                     30.21
     &         DEP2(MCGRD)          ,                                     30.21
     &         HSACC2(MCGRD)        ,                                     30.21
     &         SACC2(MCGRD)                                               30.21
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'INSAC')
!
      IF ( LMXF ) THEN                                                    40.41 40.30
         IX1 = 1
      ELSE
         IX1 = 1+IHALOX
      END IF
      IF ( LMXL ) THEN
         IX2 = MXC
      ELSE
         IX2 = MXC-IHALOX
      END IF
      IF ( LMYF ) THEN
         IY1 = I1MYC
      ELSE
         IY1 = 1+IHALOY
      END IF
      IF ( LMYL ) THEN
         IY2 = I2MYC
      ELSE
         IY2 = MYC-IHALOY
      END IF
!
      DO 30 IX = IX1, IX2
      DO 20 IY = IY1, IY2
!
!       *** Compute the mean ENERGY DENSITY frequency SACC2  ***
!       *** and the wavenumber HSACC2 average over the full  ***
!       *** spectrum per gridpoint                           ***
!
        IND = KGRPNT(IX,IY)
        IF (DEP2(IND) .GT. DEPMIN ) THEN                                  25/MAR
          SME_T  = 0.
          SME_B  = 0.
          DO 18 IS = 1, MSC
            DO 17 ID = 1, MDC
              ACS2 = SPCSIG(IS)**2 * AC2(ID,IS,IND)                       30.72
              ACS3 = SPCSIG(IS) * ACS2                                    30.72
              SME_B = SME_B + ACS2
              SME_T = SME_T + ACS3
 17         CONTINUE
 18       CONTINUE
          SME_B = SME_B * FRINTF * DDIR                                   30.50
          SME_T = SME_T * FRINTF * DDIR                                   30.50
!
!         *** mean frequency and significant wave height ***
!         *** per gridpoint                               ***
!
          IF ( SME_B .LE. 0. ) THEN
            SACC2(IND)  = 1.E-20
            HSACC2(IND) = 1.E-20
          ELSE
            SACC2(IND)  = MAX ( 1.E-20 , (SME_T / SME_B) )
            HSACC2(IND) = MAX ( 1.E-20 , (4. * SQRT(SME_B)) )
          END IF
        ELSE
          SACC2(IND)  = 0.
          HSACC2(IND) = 0.
        END IF
 20   CONTINUE
 30   CONTINUE
!
!     End of the subroutine INSAC
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE ACTION (IDCMIN     ,IDCMAX     ,SPCSIG     ,             33.09
     &                   AC2        ,CAX        ,CAY        ,
     &                   CAS        ,CAD        ,IMATLA     ,
     &                   IMATDA     ,IMATUA     ,IMATRA     ,
     &                   IMAT5L     ,                                     40.22
     &                   IMAT6U     ,ISCMIN     ,ISCMAX     ,
     &                   IDDLOW     ,IDDTOP     ,ISSTOP     ,
     &                   ANYBLK     ,ANYBIN     ,
     &                   LEAKC1     ,AC1        ,                         40.00
     &                   DYNDEP     ,RDX        ,RDY        ,             30.51
     &                   SWPDIR     ,IX         ,IY         ,
     &                   KSX        ,KSY        ,
     &                   XCGRID     ,YCGRID     ,                         40.41
     &                   ITER       ,KGRPNT     ,OBREDF     ,             40.41
     &                   CAX1       ,CAY1       ,SPCDIR     ,             33.08
     &                   CGO        ,TRAC0      ,TRAC1                    40.85 33.08
     &                                                      )
!
!****************************************************************
!
      USE TIMECOMM                                                        40.41
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
!     30.74: IJsbrand Haagsma (Include version)
!     33.08: W. Erick Rogers (some S&L scheme-related changes)
!     33.09: Nico Booij and Erick Rogers
!     33.10: Nico Booij and Erick Rogers
!     40.03: Nico Booij
!     40.08: Erick Rogers
!     40.09: Annette Kieftenburg
!     40.22: John Cazes and Tim Campbell
!     40.23: Marcel Zijlema
!     40.28: Annette Kieftenburg
!     40.41: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Mar. 98: water level (WLEV2) and wave height (CHS) in comp. grid
!                     added as arguments (needed for SWTRCF)
!                     Call SWTRCF modified
!     33.08, July 98: some S&L scheme-related changes
!     33.09, Sept 99: changes re: the spherical coordinates
!     33.10, Jan. 00: changes re: the SORDUP scheme
!     40.09, May  00: Argument list SWTRCF modified
!     40.03, Apr. 00: integers LINK1 and LINK2 replaced by array LINK(1:MICMAX)
!     40.22, Sep. 01: Removed WAREA array.                                40.22
!     40.22, Sep. 01: Changed array definitions to use the parameter      40.22
!                     MICMAX instead of ICMAX.                            40.22
!     40.28, Dec. 01: Argument list SWTRCF modified                       40.28
!     40.23, Aug. 02: Print of CPU times added
!     40.23, Nov. 02: call to SWFLXD added
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.41, Aug. 04: call to SWTRCF removed because superfluous
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.85, Aug. 08: add transport for output purposes
!
!  2. Purpose
!
!     to determine the transport, refraction and the source terms
!     of the action balance equation
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!     XCGRID: Coordinates of computational grid in x-direction            30.72
!     YCGRID: Coordinates of computational grid in y-direction            30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!     IX          Counter of gridpoints in x-direction
!     IY          Counter of gridpoints in y-direction
!     ICMAX       Maximum counter for the points of the molecul
!     MXC         Maximum counter of gridppoints in x-direction
!     MYC         Maximum counter of gridppoints in y-direction
!     MSC         Maximum counter of relative frequency
!     MDC         Maximum counter of directional distribution
!     KSX         Dummy variable to get the right sign in the
!                 numerical difference scheme in X-direction
!                 depending on the sweep direction, KSX = -1 or +1
!     KSY         Dummy variable to get the right sign in the
!                 numerical difference scheme in Y-direction
!                 depending on the sweep direction, KSY = -1 or +1
!     IDTOT,ISTOT Maximum range between the counters in directional
!                 space and frequency space respectively
!     IDDLOW      Minimum counter per sweep taken over all
!                 frequencies
!     IDDTOP      Maximum counter per sweep taken over all
!                 frequencies
!     ISSTOP      Maximum counter per sweep taken over all
!                 frequencies
!
!
!     REALS:
!     ---------
!
!     DX,DY       Step size in x-direction and y-direction
!     DDX         Same as DX but with correct sign depending of the
!                 direction of the sweep
!     DDY         Same as DY but with correct sign depending of the
!                 direction of the sweep
!
!     one and more dimensional arrays:
!     ---------------------------------
!
!     AC2       4D    Action density as function of D,S,X,Y at time T
!     CAD       3D    Wave transport velocity in spectral direction as
!                     function of (ID,IS,IC)
!     CAS       3D    Wave transport velocity in frequency-direction as
!                     as function of (ID,IS,IC)
!     CAX       3D    Wave transport velocity in X-dirction as function of
!                     (ID,IS,IC)
!     CAY       3D    Wave transport velocity in Y-dirction as function of
!                     (ID,IS,IC)
!     IMATDA    2D    Coefficients of main diagonal of matrix
!     IMATLA    2D    Coefficients of lower diagonal of matrix
!     IMATUA    2D    Coefficients of upper diagonal of matrix
!     IMATRA    2D    Coefficients of right hand side
!     IMAT5L    2D    coefficient of lower diagonal in presence of
!                     a current (see routine SWSIP)
!     IMAT6U    2D    coefficient of upper diagonal in presence of
!                     a current (see routine SWSIP)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     STRSXY
!     SORDUP
!     SANDL
!     STRSSI
!     STRSSB
!     STRSD
!     SWFLXD
!
!  9. Subroutines calling
!
!     SWOMPU
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
!     -------------------------------------------------------------------
!     *** transport in geographical space ***
!     Call STRSX to compute the propagation terms in X-direction
!     Call STRSY to compute the propagation terms in Y-direction
!     -------------------------------------------------------------------
!     *** transport in frequency space ***
!     If implicit scheme in frequency space then
!     ---
!       Call STRSSI  to compute the propagation terms in S-direction
!     ---
!     else if explicit scheme in frequency space and the energy near
!          the blocking point is removed from the spectrum then
!     ---
!       Call STRSSB to compute the propagation terms in S-direction
!     ---
!     endif
!     -------------------------------------------------------------------
!     IF no flux-limiting DO
!       Call STRSD to compute the propagation terms in directional domain
!     ELSE IF flux-limiting DO
!       Call SWFLXD
!     ---------------------------------------------------------
!     End of subroutine ACTION
!     ---------------------------------------------------------
!
!  13. Source text
!
      INTEGER  IDDLOW  ,IDDTOP  ,ISSTOP  ,                                33.09 NB!
     &                  SWPDIR  ,ITER                                     33.09 NB!
!
      LOGICAL           DYNDEP                                            33.09
!
      INTEGER :: IDCMIN(MSC), IDCMAX(MSC)
      INTEGER :: ISCMIN(MDC), ISCMAX(MDC)
      INTEGER :: KGRPNT(MXC,MYC)
!
      REAL  :: AC2(MDC,MSC,MCGRD)  ,AC1(MDC,MSC,MCGRD)
      REAL  :: CAX(MDC,MSC,MICMAX)  ,CAY(MDC,MSC,MICMAX)                  40.22
      REAL  :: CAX1(MDC,MSC,MICMAX) ,CAY1(MDC,MSC,MICMAX)                 33.08 40.22
      REAL  :: CGO(MSC,MICMAX)                                            33.08 40.22
      REAL  :: CAS(MDC,MSC,MICMAX) ,CAD(MDC,MSC,MICMAX)                   40.22
      REAL  :: IMATLA(MDC,MSC)     ,IMATDA(MDC,MSC)     ,
     &         IMATUA(MDC,MSC)     ,IMATRA(MDC,MSC)     ,
     &         IMAT5L(MDC,MSC)     ,IMAT6U(MDC,MSC)     ,
     &         LEAKC1(MDC,MSC)
      REAL  :: RDX(10)             ,RDY(10)             ,                 33.09 40.08
     &         OBREDF(MDC,MSC,2)   ,
     &         SPCDIR(MDC,6)                                              33.08
      REAL  :: TRAC0 (1:MDC,1:MSC,1:MTRNP)                                40.85
      REAL  :: TRAC1 (1:MDC,1:MSC,1:MTRNP)                                40.85
!
      LOGICAL  ANYBLK(MDC,MSC)     ,ANYBIN(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'ACTION')
!
!     *** set the coefficients in the arrays 0 ***
!
      DO IS = 1, MSC                                                      17/JAN
        DO ID = 1, MDC
          IMATLA(ID,IS) = 0.
          IMATUA(ID,IS) = 0.
          IMAT5L(ID,IS) = 0.                                              40.23
          IMAT6U(ID,IS) = 0.                                              40.23
        ENDDO
      ENDDO
!
!     set leak coefficient at 0
!
      DO ISC = 1, MSC
        DO IDC = 1, MDC
          LEAKC1(IDC,ISC) = 0.
        ENDDO
      ENDDO
!
!     *** set all transport coeff at 0 ***
!
      TRAC0(1:MDC,1:MSC,1:MTRNP) = 0.                                     40.85
      TRAC1(1:MDC,1:MSC,1:MTRNP) = 0.                                     40.85
!
!TIMG      CALL SWTSTA(140)                                                    40.23
!
!     *** Call propagation module in X-Y space  ***
!
!     --- depending on PROPSL, call STRSXY or other scheme                33.10
      IF (PROPSL.EQ.3) THEN    ! use S&L scheme
          CALL SANDL(ISSTOP   ,IDCMIN   ,IDCMAX   ,CGO     ,CAX    ,      33.09
     &               CAY      ,AC2      ,AC1      ,IMATRA  ,IMATDA ,      33.09
     &               RDX      ,RDY      ,CAX1     ,CAY1    ,SPCDIR ,      33.09
     &               TRAC0    ,TRAC1    )                                 40.85
      ELSE IF (PROPSL.EQ.2) THEN ! use SORDUP scheme                      33.10
         CALL SORDUP(ISSTOP   ,IDCMIN   ,IDCMAX   ,CAX      ,             33.10
     &               CAY      ,AC2      ,IMATRA   ,IMATDA   ,             33.10
     &               RDX      ,RDY      ,TRAC0    ,TRAC1    )             40.85 33.10
      ELSE                     ! use BSBT scheme
         CALL STRSXY(ISSTOP   ,IDCMIN   ,IDCMAX   ,CAX      ,             40.00
     &               CAY      ,AC2      ,AC1      ,IMATRA   ,IMATDA   ,
     &                         RDX      ,RDY      ,                       40.00
     &               OBREDF   ,TRAC0    ,TRAC1    )                       40.85 40.00
!
      END IF                                                              33.08
!TIMG      CALL SWTSTO(140)                                                    40.23
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 120 ) THEN
        WRITE(PRINTF,111) KCGRD(1)
 111    FORMAT(' ACTION POINT  :',I5)
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' matrix coefficients in action after strs(x-y)'
        WRITE(PRINTF,*)
        WRITE(PRINTF,*)
     & '   IS   ID    IMAT5L       IMATDA       IMAT6U    IMATRA    CAS'
        DO IDDUM = IDDLOW, IDDTOP
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IF (ISCMIN(ID).GT.0) THEN
            DO IS = ISCMIN(ID), ISCMAX(ID)
              WRITE(PRINTF,2101) IS, ID, IMAT5L(ID,IS),IMATDA(ID,IS),
     &                         IMAT6U(ID,IS),IMATRA(ID,IS),CAS(ID,IS,1)
2101          FORMAT(1X,2I4,4X,4E12.4,E10.2)
            ENDDO
          ENDIF
        ENDDO
      END IF
!
!TIMG      CALL SWTSTA(141)                                                    40.23
      IF ( (DYNDEP .OR. ICUR.EQ.1) .AND. ITFRE.NE.0 ) THEN                40.41 40.00
!
!       *** call propagation module in S-direction ***
!
        IF ( INT(PNUMS(8)) .EQ. 1 ) THEN
!
!         *** use implicit scheme for the integration in frequency ***
!         *** space (no a priori assumptions)                      ***
!
          CALL STRSSI (SPCSIG  ,
     &                 CAS     ,IMAT5L  ,IMATDA  ,IMAT6U  ,ANYBIN  ,
     &                 IMATRA  ,AC2     ,ISCMIN  ,ISCMAX  ,IDDLOW  ,
     &                 IDDTOP  ,TRAC0   ,TRAC1                     )      40.85 40.41 30.21
!
        ELSE IF ( INT(PNUMS(8)) .EQ. 2 ) THEN
!
!         *** Explicit numerical scheme in frequency space    ***
!         *** based on flux transport of action across        ***
!         *** boundaries. Energy is removed from the spectrum ***
!         *** based on a CFL criterion                        ***
!
          CALL STRSSB (IDDLOW  ,IDDTOP  ,
     &                 IDCMIN  ,IDCMAX  ,ISSTOP  ,CAX     ,CAY     ,
     &                 CAS     ,AC2     ,SPCSIG  ,IMATRA  ,
     &                 ANYBLK  ,RDX     ,RDY     ,TRAC0            )      41.07 40.41 30.21
!
        END IF
      END IF
!TIMG      CALL SWTSTO(141)                                                    40.23
!
!     *** call propagation module in D-direction ***
!
!TIMG      CALL SWTSTA(142)                                                    40.23
      IF ( IREFR.NE.0 ) THEN                                              40.41
         IF ( PROPFL.EQ.0 ) THEN                                          40.23
            CALL STRSD (DDIR    ,IDCMIN  ,
     &                  IDCMAX  ,CAD     ,IMATLA  ,IMATDA  ,IMATUA  ,
     &                  IMATRA  ,AC2     ,ISSTOP  ,
     &                  ANYBIN  ,LEAKC1  ,TRAC0   ,TRAC1            )     40.85 40.41 30.21
         ELSE IF ( PROPFL.EQ.1 ) THEN                                     40.23
            CALL SWFLXD (CAD, IMATLA, IMATDA, IMATUA, IMATRA,             40.23
     &                   AC2, DDIR, ANYBIN, LEAKC1, IDCMIN,               40.23
     &                   IDCMAX, ISSTOP)                                  40.23
         END IF                                                           40.23
      END IF                                                              40.41
!TIMG      CALL SWTSTO(142)                                                    40.23
!
!     *** test; remove on vector computer ***
!
      IF ( TESTFL .AND. ITEST .GE. 70 ) THEN
        WRITE(PRINTF,*) ' *** Values at end of subroutine action ***'
        WRITE (PRINTF,6120) KCGRD(1), MCGRD, MSC, MDC
 6120   FORMAT (' ACTION: POINT MCGRD MSC MDC   : ',4I5)
        WRITE (PRINTF,6220) IDDLOW,IDDTOP,ISSTOP, ICMAX
 6220   FORMAT (' ACTION: IDLW IDTP ISTOP ICMAX   : ',4I4)
        WRITE (PRINTF,6020) KCGRD(1), KCGRD(2), KCGRD(3)
 6020   FORMAT (' ACTION: KCGRD(1), KCGRD(2), KCGRD(3)     : ',3I4)
        WRITE (PRINTF,6022) RDX(1), RDX(2), RDY(1), RDY(2)
 6022   FORMAT (' ACTION:RDX(1) RDX(2) RDY(1) RDY(2) : ',4E12.4)
        IF (ITEST.GE.210) THEN                                            40.00
          DO IS = 1, MSC
            WRITE(PRINTF,6420) SPCSIG(IS),IDCMIN(IS),IDCMAX(IS)           30.72
 6420       FORMAT (' ACTION: SPCSIG IDCMIN IDCMAX     : ',F8.4,2I6)      30.72
          ENDDO
          DO IS = 1, MSC
            WRITE(PRINTF,*) 'IS ',IS
            DO ID = 1, MDC
              WRITE(PRINTF,6320) ID, CAX(ID,IS,1), CAY(ID,IS,1),
     &                   CAS(ID,IS,1), CAD(ID,IS,1), AC2(ID,IS,KCGRD(1))
 6320         FORMAT(' ACTION: ID CAX CAY CAS CAD AC2:',I3,5E12.4)
            ENDDO
          ENDDO
        ENDIF                                                             40.00
        WRITE(PRINTF,*) ' *** end of subr ACTION *** '                    40.00
      END IF
!     End of subroutine ACTION
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SINTGRL(SPCDIR  ,KWAVE   ,AC2     ,                      40.02
     &                   DEP2    ,QB_LOC  ,URSELL  ,                      40.02
     &                   RDX     ,RDY     ,                               40.02
     &                   AC2TOT  ,ETOT    ,                               40.02
     &                   ABRBOT  ,UBOT    ,HS      ,QB      ,             40.02
     &                   HM      ,KMESPC  ,SMEBRK  ,                      40.02
     &                   TMBOT   ,SWPDIR           )                      40.51 40.16
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
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
!     40.08: Erick Rogers
!     40.13: Nico Booij
!     40.16: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!     40.41: Marcel Zijlema
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.02, Jan. 00: New, based on the old SDISPA subroutine
!     40.02, Oct. 00: KWAVE removed in call BRKPAR
!     40.13, Aug. 01: reduction of spectrum not for Mode Noupdate
!     40.22, Sep. 01: Changed array definitions to use the parameter      40.22
!                     MICMAX instead of ICMAX.                            40.22
!     40.22, Sep. 01: Changed allocated arrays to static arrays to fix    40.22
!                     OpenMP problems with arrays allocated in parallel   40.22
!                     regions.                                            40.22
!     40.22, Oct. 01: PSURF(2) is no longer used as a variable
!     40.16, Dec. 01: Implemented limiter switches
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: near bottom wave period added
!
!  2. Purpose
!
!     To compute several integrals used in SWAN and some general parameters
!
!  3. Method
!
!     The total energy ETOT is calculate as the following integral
!
!     ETOT = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!     To avoid too high dissipation by breaking, ETOT is maximised by a maximum
!     total energy EMAX based on the maximum wave heigth HM:
!
!     HM    = PSURF(2) depth
!
!                      2
!     EMAX  = 0.25 * HM
!
!     When EMAX > ETOT, then the action density AC2 is reduced by EMAX/ETOT
!
!     In the physicaly unrealistic case that ETOT <= 0, the integrals and other parameters
!     get values that represent a steady sea-state with wind close to zero.
!
!     The following integrals are calculated:
!
!                                                       2
!     AB2   = Integrate [ (AC2(theta,sigma) sigma / Sinh [ K(sigma) depth ]) dsigma dtheta ]
!
!     ACTOT = Integrate [ AC2(theta,sigma) dsigma dtheta ]
!
!     EDRKTOT=Integrate [ (AC2(theta,sigma) sigma / Sqrt [ K(sigma) ]) dsigma dtheta ]
!
!     EKTOT = Integrate [ AC2(theta,sigma) K(sigma) sigma dsigma dtheta ]
!
!                                               2
!     ETOT1 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                               3
!     ETOT2 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                               5
!     ETOT4 = Integrate [ AC2(theta,sigma) sigma dsigma dtheta ]
!
!                                                3      2
!     UB2   = Integrate [ (AC2(theta,sigma) sigma / Sinh [ K(sigma) depth ]) dsigma dtheta ]
!
!     For reasons of ??, in the calculation of UB2, AB2, ETOTM2, ETOTM4, the high frequency
!     tail is ignored.
!
!     Based on these integrals the following parameters are calculated:
!
!     ABRBOT  = Sqrt [ 2 AB2 ]
!     HS      = 4 Sqrt [ ETOT ]
!                               2
!     KM_WAM  = (ETOT / EDRKTOT)
!     QB      : computed in the subroutine FRABRE
!     SIGM01  = ETOT1 / ETOT
!     SIGM_10 = ETOT  / ACTOT
!     UBOT    = Sqrt [ UB2 ]             NOTE: THIS IS THE ROOT MEAN SQUARE OF THE ORBITAL MOTION NEAR THE BOTTOM!!!
!     TMBOT   = PI Sqrt [ 2 AB2 / UB2 ]
!
!  4. Argument variables
!
!     ABRBOT: Near bottom excursion
!     AC2   : Action density as function of ID, IS, IX and IY
!     AC2TOT : Total action density per gridpoint
!     DEP2  : Water depth
!     ETOT  : Total wave energy density
!     HM    : Maximum wave height
!     HS    : Significant wave height
!     KMESPC: Mean average wavenumber according to the WAM-formulation
!     KWAVE : Wavenumber function of frequency and IC??
!     QB    : Fraction of breaking waves
!     QB_LOC: Fraction of breaking waves at current grid-point
!     SMEBRK: Mean frequency according to first order moment
!     SWPDIR: Number of current sweep direction                           40.16
!     UBOT  : Near bottom velocity as function of IX and IY
!     URSELL: Ursell number as function of IX and IY
!     TMBOT : near bottom wave period                                     40.51

      INTEGER, INTENT(IN)  :: SWPDIR                                      40.16

      REAL, INTENT(IN)     :: DEP2(MCGRD)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL, INTENT(IN)     :: KWAVE(MSC,MICMAX)                           40.22
      REAL, INTENT(IN)     :: RDX(10), RDY(10)                            40.08
      REAL, INTENT(IN)     :: SPCDIR(MDC,6)
!
      REAL, INTENT(IN OUT) :: AC2(MDC,MSC,MCGRD)
      REAL, INTENT(IN OUT) :: QB(MCGRD)
      REAL, INTENT(IN OUT) :: UBOT(MCGRD)
      REAL, INTENT(IN OUT) :: URSELL(MCGRD)
      REAL, INTENT(IN OUT) :: TMBOT(MCGRD)                                40.51
!
      REAL, INTENT(OUT)    :: ABRBOT, ETOT, HM, HS, QB_LOC
      REAL, INTENT(OUT)    :: AC2TOT, KMESPC, SMEBRK
!
!  6. Local variables
!
!     IENT  : Number of entries in this subroutine
!     IS    : Counter for the relative frequency band
!
      INTEGER, SAVE :: IENT = 0
!
!     AB2                : Sum of E_DSHKD2_DS_DD
!     ACTOT_DSIG         : Integration term for calculating ACTOT
!     EDRKTOT_SWELL      : Swell-part of EDRKTOT
!     EMAX               : Maximum energy according calculated HM
!     ETOT_DRK_DSIG      : Integration term for calculating EDRKTOT
!     ETOT_K_DSIG        : Integration term for calculating EKTOT
!     ETOT_DSHKD2_DSIG   : Integration term for calculating AB2
!     ETOT_DSIG          : Integration term for calculating ETOT
!     ETOT_SIG_DSIG      : Integration term for calculating ETOT1
!     ETOT_SIG2_DSIG     : Integration term for calculating ETOT2
!     ETOT_SIG2_DSHKD2_DSIG: Integration term for calculating UB2
!     ETOT_SIG4_DSIG     : Integration term for calculating ETOT4
!     FRINT_X_DDIR       : FRINTF * DDIR
!     SINH_K_X_DEP_2     : SINH(KWAVE*DEP2)**2
!     UB2                : Sum of ETOT_SIG2_DSHKD2_DSIG
!     WH                 : fraction of breaking waves
!
      REAL              :: AB2, EMAX, UB2
      REAL              :: WH
      REAL              :: FRINTF_X_DDIR
      REAL              :: BRCOEF   ! variable breaking coefficient (calc. in BRKPAR)  40.22
!
      REAL              :: ETOT_DSIG(MSC)                                 40.22
      REAL              :: ACTOT_DSIG(MSC), ETOT_DSHKD2_DSIG(MSC)         40.22
      REAL              :: ETOT_DRK_DSIG(MSC), ETOT_SIG2_DSIG(MSC)        40.22
      REAL              :: ETOT_K_DSIG(MSC), ETOT_SIG_DSIG(MSC)           40.22
      REAL              :: ETOT_SIG2_DSHKD2_DSIG(MSC)                     40.22
      REAL              :: ETOT_SIG4_DSIG(MSC)                            40.22
      REAL              :: SINH_K_X_DEP_2(MSC)                            40.22
!
!     9. STRUCTURE
!
!   ----------------------------------------------------------
!   Determine ETOT
!   If energy level is too large compared with depth
!   Then reduce action densities
!   ------------------------------------------------------
!   For all spectral frequencies do
!       determine wavenumber K and K*depth
!       For every spectral direction do
!           add AC2 to sum of action densities
!       --------------------------------------------------
!       add contributions to various moments of energy density
!   ---------------------------------------------------------
!   add tail contributions
!   determine average frequency and wavenumber
!   ----------------------------------------------------------
!   If B&J surf breaking is used
!   Then call FRABRE to compute fraction of breaking waves
!   ----------------------------------------------------------
!   determine orbital motion near the bottom
!   ----------------------------------------------------------
!
! 13. Source text:
!
      IF (LTRACE) CALL STRACE (IENT,'SINTGRL')
!
!     Initialisation
!
      KM_WAM  = 10.
      KM01    = 10.
!
      SIGM01  = 10.
      SIGM_10 = 10.
!
      HS      = 0.
      HM      = 0.1
!
      QB(KCGRD(1))   = 0.
      ABRBOT         = 0.001
      UBOT(KCGRD(1)) = 0.
      TMBOT(KCGRD(1))= 0.
!
!     Calculate total spectral energy:
!
!
      FRINTF_X_DDIR = FRINTF * DDIR
      ETOT_DSIG(:)  = SUM(AC2(:,:,KCGRD(1)),DIM=1) * SIGPOW(:,2) *
     &                FRINTF_X_DDIR
      ETOT          = SUM(ETOT_DSIG)
!
!     *** add high frequency tail ***
!
      ETOT = ETOT + ETOT_DSIG(MSC) * PWTAIL(6) / FRINTF
!
      IF (ISURF .EQ. 1 .OR. ISURF.EQ.4                                    41.03
     &                                ) THEN
        HM   = PSURF(2) * DEP2(KCGRD(1))                                  40.22
      ELSEIF (ISURF .EQ. 2 .OR. ISURF.EQ.3) THEN                          41.03
!        Calulate the correct breaking coefficient BRCOEF
         CALL BRKPAR (BRCOEF  ,SPCDIR(1,2), SPCDIR(1,3), AC2     ,
     &                SIGPOW(:,1), DEP2, RDX, RDY, KWAVE         )        41.03 40.22
        HM   = BRCOEF * DEP2(KCGRD(1))                                    40.22
      ELSE
!       breaking disabled, assign very high value to HM
        HM   = 100.                                                       40.22
      ENDIF
!
      EMAX = 0.25 * HM**2
!
!     --- reduce action density if necessary
!
      IF (ACUPDA .AND. ETOT .GT. EMAX .AND. ISURF .GE. 1) THEN            40.13
         AC2(:,:,KCGRD(1)) = MAX(0.,(EMAX/ETOT)*AC2(:,:,KCGRD(1)))
!
         IF (TESTFL.AND.ITEST.GE.80)
     &      WRITE (PRTEST,7) DEP2(KCGRD(1)), EMAX, ETOT
   7     FORMAT (' energy is reduced in SINTGRL', 4(1x, e12.4))
!
!        Correct value for ETOT
!
         ETOT = EMAX
!
      ENDIF
!
      IF ( ETOT .GT. 0. ) THEN
!
!     Calculate all other integrals
!
!
        SINH_K_X_DEP_2(:)   = SINH(MIN(30.,
     &                             KWAVE(:,1)*DEP2(KCGRD(1)))
     &                            )**2
        ACTOT_DSIG(:)       = SUM(AC2(:,:,KCGRD(1)),DIM=1) *
     &                        SIGPOW(:,1) * FRINTF_X_DDIR
        ETOT_DSIG(:)        = ACTOT_DSIG(:) * SIGPOW(:,1)
        ETOT_SIG_DSIG(:)    = ACTOT_DSIG(:) * SIGPOW(:,2)
        ETOT_SIG2_DSIG(:)   = ACTOT_DSIG(:) * SIGPOW(:,3)
        ETOT_SIG4_DSIG(:)   = ACTOT_DSIG(:) * SIGPOW(:,5)
        ETOT_DRK_DSIG(:)    = ETOT_DSIG(:) / SQRT(KWAVE(:,1))
        ETOT_K_DSIG(:)      = ETOT_DSIG(:) * KWAVE(:,1)
        ETOT_DSHKD2_DSIG(:) = ETOT_DSIG(:) / SINH_K_X_DEP_2(:)
        ETOT_SIG2_DSHKD2_DSIG(:) = ETOT_DSHKD2_DSIG(:) *
     &                             SIGPOW(:,2)
!
        ACTOT         = SUM(ACTOT_DSIG)
        ETOT1         = SUM(ETOT_SIG_DSIG)
        ETOT2         = SUM(ETOT_SIG2_DSIG)
        ETOT4         = SUM(ETOT_SIG4_DSIG)
        EDRKTOT       = SUM(ETOT_DRK_DSIG)
        EKTOT         = SUM(ETOT_K_DSIG)
        UB2           = SUM(ETOT_SIG2_DSHKD2_DSIG)
        AB2           = SUM(ETOT_DSHKD2_DSIG)
!
!     Add high frequency tails
!
        ACTOT       = ACTOT + PWTAIL(5) * ACTOT_DSIG(MSC) / FRINTF
        ETOT1       = ETOT1 + PWTAIL(7) * ETOT_DSIG(MSC) *
     &                SIGPOW(MSC,1) / FRINTF
        EDRKTOT     = EDRKTOT + PWTAIL(5) * ETOT_DSIG(MSC) /
     &                (SQRT(KWAVE(MSC,1)) * FRINTF)
        EKTOT       = EKTOT + PWTAIL(8) * ETOT_DSIG(MSC) *
     &                KWAVE(MSC,1) / FRINTF
!
!     Calculate the mean frequencies SIGM01 and SIGM_10,
!     mean wavenumbers KM_WAM, KM01 and significant waveheight HS
!
        IF (ETOT1 .GT. 0.) SIGM01  = ETOT1 / ETOT
        IF (EKTOT .GT. 0.) KM01    = EKTOT / ETOT
        IF (ACTOT .GT. 0.) SIGM_10 = ETOT / ACTOT
        IF (EDRKTOT .GT. 0. ) THEN
          KM_WAM  = ( ETOT / EDRKTOT )**2.
        ENDIF
        IF ( ETOT .GT. 1.E-20 ) THEN
          HS       = 4. * SQRT (ETOT)
        END IF
!
!       --- calculate Qb when Battjes/Janssen breaking is activated
!
        IF ( ISURF.GT.0 .AND. ISURF.LT.4 ) THEN
           CALL FRABRE (HM, ETOT, QB(KCGRD(1)))
!
!       --- calculate Qb when Thornton/Guza breaking is activated         41.03
!
        ELSEIF (ISURF.EQ.4) THEN
           WH = (2.*SQRT(2.*ETOT)/HM)**PSURF(4)
           WH = MIN(1.,WH)
           QB(KCGRD(1)) = WH
        ENDIF
!
!     Calculate the orbital velocity UBOT, orbital excursion ABRBOT and
!     near bottom wave period TMBOT                                       40.51
!
        IF ( UB2 .GT. 0.) UBOT(KCGRD(1)) = SQRT ( UB2 )
        IF ( AB2 .GT. 0.) ABRBOT = SQRT (2. *  AB2)
        IF ( UB2.GT.0. .AND. AB2.GT.0. )                                  40.51
     &                             TMBOT(KCGRD(1)) = PI*SQRT(2.*AB2/UB2)  40.51
!
      ENDIF
!
!     *** calculate Ursell number ***
!     *** update only for first encounter in a sweep                      40.16
!
      IF ( ITRIAD.GT.0                                                    40.41
     &                 ) THEN                                             40.41
         IF (( SWPDIR .EQ. 1) .OR.                                        40.16
     &       ( SWPDIR .EQ. 2 .AND. IXCGRD(1) .EQ. 1) .OR.                 40.16
     &       ( SWPDIR .EQ. 3 .AND. IYCGRD(1) .EQ. 1) .OR.                 40.16
     &       ( SWPDIR .EQ. 4 .AND.                                        40.16
     &           (IXCGRD(1).EQ.MXC .AND. IYCGRD(1).EQ.1) )) THEN          40.16
           URSELL(KCGRD(1)) = (GRAV*HS) /
     &                     (2.*SQRT(2.)*SIGM01**2*DEP2(KCGRD(1))**2)
         ENDIF                                                            40.16
      ELSE                                                                40.41
        URSELL(KCGRD(1)) = 0.                                             40.41
      END IF                                                              40.41
!
      QB_LOC = QB(KCGRD(1))
!
!     *** test output ***
!
      IF (TESTFL .AND. ITEST.GE.60) THEN
         WRITE(PRTEST, 901) ETOT, HS, SIGM_10, KM_WAM, ABRBOT
 901     FORMAT (' SINTGRL: ETOT Hs Sigma K Aorb', 5(1X, E11.4))
      END IF
!
!     Set variables used outside the whitecapping scope
!
      AC2TOT = ACTOT
      KMESPC = KM_WAM
      SMEBRK = SIGM01
!
      RETURN
!
      END SUBROUTINE SINTGRL
!
!****************************************************************
!
      SUBROUTINE SOLPRE (AC2         ,AC2OLD      ,                       40.00
     &                   IMATRA      ,IMATLA      ,
     &                   IMATDA      ,IMATUA      ,
     &                   IMAT5L      ,IMAT6U      ,
     &                   IDCMIN      ,IDCMAX      ,
     &                   ANYBIN      ,
     &                   IDTOT       ,ISTOT       ,
     &                   IDDLOW      ,IDDTOP      ,
     &                   ISSTOP      ,
     &                   SPCSIG                   )                       40.41 40.23
!
!****************************************************************
!
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.00: Nico Booij
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New subroutine common tasks before solution of linear
!                     system (software moved from SOLBAND, SOLMAT and SOLMT1)
!     40.23, Aug. 02: implementation of under-relaxation technique
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.41, Aug. 04: code optimized
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copy local spectrum to array AC2OLD before solving system,
!     apply under-relaxation approach in active bins, fill matrix
!     arrays for non-active bins and write test output
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!        AC2OLD    2D    Values of action density at previous iteration
!        IMATDA    2D    Coefficients of diagonal of matrix
!        IMATLA    2D    Coefficients of lower diagonal of matrix
!        IMATUA    2D    Coefficients of upper diagonal of matrix
!        IMATRA    2D    Coefficients of right hand side of matrix
!        IMAT5L    2D    Coefficients for implicit calculation in
!                        frequency space (lower diagonal)
!        IMAT6U    2D    Coefficients for implicit calculation in
!                        frequency space (upper diagonal)
!        SPCSIG    1D    Relative frequencies in sigma-space
!
      REAL     AC2(MDC,MSC,MCGRD)           ,
     &         IMATRA(MDC,MSC)              ,
     &         IMATLA(MDC,MSC)              ,
     &         IMATDA(MDC,MSC)              ,
     &         IMATUA(MDC,MSC)              ,
     &         IMAT5L(MDC,MSC)              ,
     &         IMAT6U(MDC,MSC)              ,
     &         AC2OLD(MDC,MSC)              ,
     &         SPCSIG(MSC)
!
!        IDCMIN    1D    Integer array containing minimum counter
!        IDCMAX    1D    Integer array containing maximum counter
!
      INTEGER  IDCMIN(MSC)                  ,
     &         IDCMAX(MSC)
!
!        ANYBIN    2D    Logical array. if a certain bin is enclosed
!                        in a sweep then ANYBIN is TRUE . array is
!                        used to determine whether some coefficients
!                        in the array have to be changed
!
      LOGICAL  ANYBIN(MDC,MSC)
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
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  IS      ,ID      ,IDDUM   ,
     &         IDDLOW  ,
     &         IDDTOP  ,IDTOT   ,ISTOT   ,ISSTOP
      REAL     ALFA                                                       40.23
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SOLPRE')
!
!     --- apply under-relaxation approach, if requested
!
      ALFA = PNUMS(30)                                                    40.41
      IF (ALFA.GT.0.) THEN                                                40.41
         DO IS = 1, ISSTOP
            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
               ID = MOD(IDDUM-1 + MDC, MDC) + 1
               IMATDA(ID,IS) = IMATDA(ID,IS) + ALFA*SPCSIG(IS)            40.23
               IMATRA(ID,IS) = IMATRA(ID,IS) + ALFA*SPCSIG(IS)*           40.23
     &                                         AC2(ID,IS,KCGRD(1))        40.23
            END DO
         END DO
      END IF
!
!     --- when ambient currents are involved or when the spectral space
!         is not a circular one (use SECTOR instead of CIRCLE in command
!         CGRID), some bins do not fall within the current sweep (in
!         particular, when SECTOR = 0 or 4 i.c. ICUR=1 or SECTOR = 4 i.c.
!         FULCIR=.FALSE., see routine SWPSEL for meaning of SECTOR). For
!         such bins, the corresponding rows in the matrix are reset such
!         that the solution AC2 does not change: the main diagonal is set
!         to 1, the off-diagonals are set to 0 and the righ-hand side is
!         set to AC2.
!
      IF ( ICUR.EQ.1 .OR. .NOT.FULCIR ) THEN                              40.41
         DO IS = 1, MSC
           DO ID = 1, MDC
             IF ( .NOT. ANYBIN(ID,IS) ) THEN
               IMATLA(ID,IS) = 0.
               IMATDA(ID,IS) = 1.
               IMATUA(ID,IS) = 0.
               IMATRA(ID,IS) = AC2(ID,IS,KCGRD(1))                        30.21
               IMAT5L(ID,IS) = 0.
               IMAT6U(ID,IS) = 0.
             END IF
           ENDDO
         ENDDO
      END IF
!
!     *** the action density is stored in an auxiliary array AC2OLD ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          AC2OLD(ID,IS) = AC2(ID,IS,KCGRD(1))
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 70 ) THEN
        WRITE (PRINTF,6120) IXCGRD(1)+MXF-2, IYCGRD(1)+MYF-2              40.30
 6120   FORMAT (' SOLPRE: Matrix values for point:', 2I5)
        WRITE (PRINTF,6121)
 6121   FORMAT ('  bin     diagonal     r.h.s.        ID-1        ID+1',
     &          '         IS-1         IS+1')
!
        DO IS = 1, ISSTOP
          ID_MIN = IDCMIN(IS)
          ID_MAX = IDCMAX(IS)
          DO IDDUM = ID_MIN, ID_MAX
            ID = MOD(IDDUM-1 + MDC, MDC) + 1
            IF ( DYNDEP .OR. ICUR .EQ. 1 ) THEN                           40.00
              WRITE(PRINTF,6620) ID, IS, IMATDA(ID,IS), IMATRA(ID,IS),
     &        IMATLA(ID,IS), IMATUA(ID,IS), IMAT5L(ID,IS), IMAT6U(ID,IS)
 6620         FORMAT(2I3,6(1X,E12.4))
            ELSE
              WRITE(PRINTF,6620) ID, IS, IMATDA(ID,IS), IMATRA(ID,IS),
     &        IMATLA(ID,IS), IMATUA(ID,IS)
            ENDIF
          ENDDO
        ENDDO
      END IF
!
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE SOLMAT (IDCMIN     ,IDCMAX     ,                         40.00
     &                   AC2        ,IMATRA     ,
     &                   IMATDA     ,IMATUA     ,
     &                   IMATLA                                           40.23
     &                                           )
!
!****************************************************************

      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
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
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, Feb. 99: swcomm3 introduced
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     SUBROUTINE to solve the linear system which is filled in the
!     subroutine ACTION. The solution give the values for the
!     wave action for every frequency and every direction.
!     The system is solved by means of the Thomas sweep algorithm
!     in the spectral direction only.
!
!  3. Method
!
!     Solver for tri-diagonal matrix:
!
!
!        / 2  3          \ /   \
!        | 1  2  3       | |   |
!        |    1  2  3    | | N | =  RHS
!        |       1  2  3 | |   |
!        \          1  2 / \   /
!
!
!     This method consists of forward and backward sweeps.
!
!  4. Argument variables
!
!
!     IX          Counter of gridpoints in x-direction
!     IY          Counter of gridpoints in y-direction
!     IS          Counter of relative frequency band
!     ID          Counter of directional distribution
!     J           Dummy counter
!     MXC         Maximum counter of gridppoints in x-direction
!     MYC         Maximum counter of gridppoints in y-direction
!     MSC         Maximum counter of relative frequency
!     MDC         Maximum counter of directional distribution
!
!     REALS:
!     ---------
!
!     SP          Dummy variable
!     TEMP        Dummy variable
!
!     one and more dimensional arrays:
!     ---------------------------------
!     AC2       4D    Action density as function of D,S,X,Y and T
!     IMATDA    2D    Coefficients of diagonal of matrix
!     IMATLA    2D    Coefficients of lower diagonal of matrix
!     IMATUA    2D    Coefficients of upper diagonal of matrix
!     IMATRA    2D    Coefficients of right hand side of matrix
!     IDCMIN    1D    Integer array containing minimum counter
!     IDCMAX    1D    Integer array containing maximum counter
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
!     SWOMPU
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
!     -------------------------------------------------------------
!     For every D-direction within the sector do
!       Eliminate the lower diagonal
!     -------------------------------------------------------------
!     For every D-direction within the sector do
!       Solve the linear equation to get the wave action for every
!       direction (ID)
!     -------------------------------------------------------------
!     Set all the values in the arrays to zero:
!     IMATRA(MDC,MSC),IMATLA(MDC,MSC),IMATDA(MDC,MSC),IMATUA(MDC,MSC)
!     ------------------------------------------------------------
!
! 13. Source text
!
!
      INTEGER  IS      ,ID      ,J       ,ID_MIN  ,ID_MAX                 40.00
!
      REAL     SP      ,TEMP                                              40.23
!
      REAL     AC2(MDC,MSC,MCGRD)           ,
     &         IMATRA(MDC,MSC)              ,
     &         IMATLA(MDC,MSC)              ,
     &         IMATDA(MDC,MSC)              ,
     &         IMATUA(MDC,MSC)
!
      INTEGER  IDCMIN(MSC)        ,
     &         IDCMAX(MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SOLMAT')
!
!     **** 17/JAN   IN MOD (  , ) ;   + MDC WAS ADDED ****
!
      DO 180 IS = 1, MSC
        ID_MIN = IDCMIN(IS)                                               3/MAR
        ID_MAX = IDCMAX(IS)                                               3/MAR
!
!       *** elimination of the lower diagonal of the first matrix ***
!
        DO 100 IDDUM = (ID_MIN+1), ID_MAX                                 20.43
          ID   = MOD(IDDUM-1+MDC, MDC) + 1                                17/JAN
          IDM1 = MOD(IDDUM-2+MDC, MDC) + 1                                17/JAN
          SP   = IMATDA(IDM1,IS)                                          20.43
          IF ( ABS(SP) .LE. 1.E-20 ) THEN
            TEMP = IMATLA(ID,IS) / SIGN( 1.E-20 , SP)
          ELSE
            TEMP = IMATLA(ID,IS) / SP
          END IF
          IMATDA(ID,IS) = IMATDA(ID,IS) - TEMP * IMATUA(IDM1,IS)          20.43
          IMATRA(ID,IS) = IMATRA(ID,IS) - TEMP * IMATRA(IDM1,IS)          20.43
 100    CONTINUE
!
!       *** solving of the linear equations for the wave action ***
!
!       *** first for ID_MAX, then for the others ***
!
        ID   = MOD(ID_MAX-1+MDC, MDC) + 1                                 17/JAN
        SP = IMATDA(ID,IS)
        IF ( ABS(SP) .LE. 1.E-20 ) THEN
          TEMP = SIGN (1.E-20 , SP)
        ELSE
          TEMP = SP
        END IF
!
!       *** wave action for ID_MAX ***
!
        AC2(ID,IS,KCGRD(1)) = IMATRA(ID,IS) / TEMP
!
        DO 150 J = 1, (ID_MAX-ID_MIN)
          ID   = MOD(ID_MAX-J-1+MDC, MDC) +1                              17/JAN
          IDP1 = MOD(ID_MAX-J+MDC, MDC) +1                                17/JAN
          SP = IMATDA(ID,IS)
          IF ( ABS(SP) .LE. 1.E-20 ) THEN
            TEMP = SIGN (1.E-20 , SP)
          ELSE
            TEMP = SP
          END IF
          AC2(ID,IS,KCGRD(1)) = ( IMATRA(ID,IS) - IMATUA(ID,IS) *
     &                        AC2(IDP1,IS,KCGRD(1)) ) / TEMP                 20.
!
 150    CONTINUE
!
        IF ( ITEST .GE. 120 .AND. TESTFL ) THEN
          WRITE(PRINTF,6009) KCGRD(1),ID_MIN,ID_MAX
 6009     FORMAT(' SOLMAT: POINT ID_MIN ID_MAX :',3I5)
          DO 169 IDDUM = ID_MIN, ID_MAX
            ID = MOD(IDDUM-1+MDC, MDC) + 1                                17/
            WRITE (PRINTF,6010) IS,ID,AC2(ID,IS,KCGRD(1))
 6010       FORMAT(' IS ID AC2()         :',2I5,2X,E12.4)
 169      CONTINUE
        END IF
!
 180  CONTINUE
!
!     *** set all the coefficients in the arrays 0 ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          IMATRA(ID,IS) = 0.
          IMATDA(ID,IS) = 0.
        ENDDO
      ENDDO
!
      IF ( TESTFL .AND. ITEST.GE. 40 ) THEN
        WRITE (PRINTF,6020) IXCGRD(1)+MXF-2, IYCGRD(1)+MYF-2
 6020   FORMAT(' SOLMAT: point :',2I5)
        WRITE(PRINTF,*)
      END IF
!
!     End of the subroutine SOLMAT
!
      RETURN
      END
!
!****************************************************************
!
      SUBROUTINE SOLMT1 (IDCMIN     ,IDCMAX     ,                         40.00
     &                   AC2        ,IMATRA     ,
     &                   IMATDA     ,IMATUA     ,
     &                   IMATLA     ,                                     40.23
     &                   ISSTOP     ,                                     40.41 40.23
     &                   ANYBLK     ,IDDLOW     ,
     &                   IDDTOP              )
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, Feb. 99: swcomm3 introduced
!     40.41, Aug. 04: array SECTOR removed and some corrections if SECTOR=0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     subroutine to solve the linear system which is filled in the
!     subroutines ACTION and SOURCE. The solution give the values
!     for the wave action for every frequency and every direction.
!     Ambient currents are involved and the propagation term in           40.41
!     frequency space is discretized explicitly.                          40.41
!     The system is solved by means of the Thomas sweep algorithm
!     in the spectral direction only.
!
!  3. Method
!
!     Solver for tridiagonal matrix with a possible coefficient
!     at the bottom left position and top right position due
!     to periodicity in theta-space:
!
!
!        / 2  3        1 \ /   \
!        | 1  2  3       | |   |
!        |    1  2  3    | | N | =  RHS
!        |       1  2  3 | |   |
!        \ 3        1  2 / \   /
!
!
!     This method consists of forward and backward sweeps.
!
!  4. Argument variables
!
!
!        IX          Counter of gridpoints in x-direction
!        IY          Counter of gridpoints in y-direction
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        J           Dummy counter
!        MXC         Maximum counter of gridppoints in x-direction
!        MYC         Maximum counter of gridppoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        REALS:
!        ---------
!
!        SP          Dummy variable
!        TEMP        Dummy variable
!
!        one and more dimensional arrays:
!        ---------------------------------
!        AC2       4D    Action density as function of D,S,X,Y and T
!        IMATDA    2D    Coefficients of diagonal of matrix
!        IMATLA    2D    Coefficients of lower diagonal of matrix
!        IMATUA    2D    Coefficients of upper diagonal of matrix
!        IMATRA    2D    Coefficients of right hand side of matrix
!        IDCMIN    1D    Integer array containing minimum counter
!        IDCMAX    1D    Integer array containing maximum counter
!        ICOLU2    1D    In presence of a current the spectral direction can
!                        be circular and closed. Matrix coefficients appear in
!                        the top right and bottom left corner of the matrix
!                        After pivoting --> coefficients are stored in ICOLU2
!                        space
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
!   -------------------------------------------------------------
!   For every D-direction within the sector do
!     Eliminate the lower diagonal
!   -------------------------------------------------------------
!   For every D-direction within the sector do
!     Solve the linear equation to get the wave action for every
!     direction (ID)
!   -------------------------------------------------------------
!   Set all the values in the arrays to zero:
!   IMATRA(MDC,MSC),IMATLA(MDC,MSC),IMATDA(MDC,MSC),IMATUA(MDC,MSC)
!   ------------------------------------------------------------
!   End of SOLMT1
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  IS     ,ID     ,J      ,IDDUM  ,IIDM   ,
     &         IIDP   ,ISSTOP ,IDDTOP ,IDDLOW                             40.00
!
      REAL     SP     ,TEMP   ,CORMAT ,TEMP1                              40.23
!
      REAL     AC2(MDC,MSC,MCGRD)           ,                             30.21
     &         IMATRA(MDC,MSC)              ,
     &         IMATLA(MDC,MSC)              ,
     &         IMATDA(MDC,MSC)              ,
     &         IMATUA(MDC,MSC)              ,
     &         ICOLU2(MDC)
!
      INTEGER  IDCMIN(MSC)                  ,
     &         IDCMAX(MSC)
!
      LOGICAL  ANYBLK(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SOLMT1')
!
!     *** since explicit scheme is used and when CFL exceeds CFL max ***
!     *** then bin should not be propagated within the current sweep ***
!
      DO IS = 1, ISSTOP
        DO IDDUM = IDDLOW, IDDTOP
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IF ( ANYBLK(ID,IS) ) THEN
            IMATLA(ID,IS) = 0.
            IMATDA(ID,IS) = 1.
            IMATUA(ID,IS) = 0.
            IMATRA(ID,IS) = 0.
          END IF
        ENDDO
      ENDDO
!
!     *** start process of elimination ***
!
      DO 180 IS = 1, MSC
!
!         *** set values in auxiliary array for last ***
!         *** column equal zero                      ***
!
          DO ID = 1, MDC
           ICOLU2(ID) = 0.
          ENDDO
!
          IF ( IDCMIN(IS).LE.IDCMAX(IS) ) THEN                            40.41
             IDLOW = IDCMIN(IS)
             IDTOP = IDCMAX(IS)
          ELSE                                                            40.41
             IDLOW = 1
             IDTOP = MDC
          END IF                                                          40.41
!
!         *** set values in coefficients in left bottom element ***
!         *** and right top element. Situation only occurs if   ***
!         *** the matrix is solved for all directions           ***
!
          IF ( IDLOW .EQ. 1  .AND.  IDTOP .EQ. MDC ) THEN
            CORMAT    = IMATUA(MDC,IS)
            ICOLU2(1) = IMATLA(1,IS)
          ELSE
            CORMAT    = 0.
            ICOLU2(1) = 0.
          END IF
!
!         *** elimination of the lower diagonal of the first matrix ***
!
          DO 100 IDDUM = (IDLOW+1 ), IDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IIDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            SP = IMATDA(IIDM,IS)
            IF ( ABS(SP) .LE. 1.E-20 ) THEN
              TEMP = IMATLA(ID,IS) / SIGN( 1.E-20 , SP)
              TEMP1 = CORMAT / SIGN(1.E-20 , SP)
            ELSE
              TEMP = IMATLA(ID,IS) / SP
              TEMP1 = CORMAT / SP
            END IF
            IMATDA(ID,IS)  = IMATDA(ID,IS)  - TEMP * IMATUA(IIDM,IS)
            IMATRA(ID,IS)  = IMATRA(ID,IS)  - TEMP * IMATRA(IIDM,IS)
            IMATRA(IDTOP,IS) = IMATRA(IDTOP,IS) - TEMP1 *
     &                          IMATRA(IIDM,IS)
            CORMAT = 0. - TEMP1 * IMATUA(IIDM,IS)
!
            IF ( IDDUM .LT. (IDTOP-1) ) THEN
              ICOLU2(ID) =  - TEMP * ICOLU2(IIDM)
            ELSE
              IMATUA(ID,IS) = IMATUA(ID,IS) - TEMP * ICOLU2(IIDM)
            END IF
            IF ( IDDUM .LT. IDTOP ) THEN
              IMATDA(IDTOP,IS) = IMATDA(IDTOP,IS) - TEMP1 *
     &                            ICOLU2(IIDM)
            ELSE
              IMATDA(IDTOP,IS) = IMATDA(IDTOP,IS) - TEMP1 *
     &                            IMATUA(IIDM,IS)
            END IF
!
 100      CONTINUE
!
!         *** solving of the linear equations for the wave action ***
!
!         *** first for IDTOP, then for the others ***
!
          SP = IMATDA(IDTOP,IS)
          IF ( ABS(SP) .LE. 1.E-20 ) THEN
            TEMP = SIGN (1.E-20 , SP)
          ELSE
            TEMP = SP
          END IF
!
!         *** wave action for IDCMAX ***
!
          AC2(IDTOP,IS,KCGRD(1)) = IMATRA(IDTOP,IS) / TEMP
!
          DO 150 J = 1, (IDTOP-IDLOW)
            IDDUM = IDTOP - J
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IIDP = MOD ( IDDUM + MDC , MDC ) + 1
            SP = IMATDA(ID,IS)
            IF ( ABS(SP) .LE. 1.E-20 ) THEN
              TEMP = SIGN (1.E-20 , SP)
            ELSE
              TEMP = SP
            END IF
            AC2(ID,IS,KCGRD(1)) = ( IMATRA(ID,IS) - IMATUA(ID,IS) *       30.21
     &                        AC2(IIDP,IS,KCGRD(1)) - ICOLU2(ID)  *
     &                        AC2(IDTOP,IS,KCGRD(1)) ) / TEMP
 150      CONTINUE
!
!         *** extended info for SOLMT1 ***
!
          IF ( ITEST .GE. 13 .AND. TESTFL ) THEN
            WRITE(PRINTF,*) 'SOLMT1'
            WRITE(PRINTF,*) ' matrix coefficients after pivoting '
            WRITE(PRINTF,*)
            WRITE(PRINTF,*)
     &   'ID IDDUM IMATLA      IMATDA      IMATUA     ICOLU2    IMATRA'
            DO 2100 IDDUM = IDLOW, IDTOP
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
              WRITE(PRINTF,2101) ID, IDDUM,IMATLA(ID,IS),IMATDA(ID,IS),
     &                         IMATUA(ID,IS),ICOLU2(ID),IMATRA(ID,IS)
2101          FORMAT(2I3,5E12.4)
2100        CONTINUE
            WRITE(PRINTF,*)
            DO 2010 IDDUM = IDLOW, IDTOP
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
              WRITE (PRINTF,6010) IS,ID,AC2(ID,IS,KCGRD(1))
 6010         FORMAT(' IS ID and resolved vector  :',2I5,2X,E12.4)
 2010       CONTINUE
            WRITE(PRINTF,*)
          END IF
 180  CONTINUE
!
!     *** set all the coefficients in the arrays 0 ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          IMATRA(ID,IS) = 0.
          IMATDA(ID,IS) = 0.
          ICOLU2(ID)    = 0.
        ENDDO
      ENDDO
!
!     End of the subroutine SOLMT1
!
      RETURN
      END
!
!****************************************************************
!
        SUBROUTINE SOURCE (ITER       ,IX         ,IY         ,
     &                     SWPDIR     ,KWAVE      ,SPCSIG     ,           30.72
     &                     ECOS       ,ESIN       ,AC2        ,
     &                     DEP2       ,IMATDA     ,IMATRA     ,
     &                     ABRBOT     ,KMESPC     ,SMESPC     ,
     &                     UBOT       ,UFRIC      ,UX2        ,
     &                     UY2        ,IDCMIN     ,IDCMAX     ,
     &                     IDDLOW     ,IDDTOP     ,IDWMIN     ,
     &                     IDWMAX     ,ISSTOP     ,PLWNDA     ,
     &                     PLWNDB     ,PLWCAP     ,PLBTFR     ,
     &                     PLWBRK     ,PLNL4S     ,PLNL4D     ,
     &                     PLVEGT     ,                                   40.55
     &                     PLTRI      ,            HS         ,           40.22
     &                     ETOT       ,QBLOC      ,THETAW     ,
     &                     HM         ,FPM        ,WIND10     ,
     &                     ETOTW      ,GROWW      ,ALIMW      ,
     &                     SMEBRK     ,SNLC1      ,
     &                     DAL1       ,DAL2       ,DAL3       ,
     &                     UE         ,SA1        ,                       40.17
     &                     SA2        ,DA1C       ,DA1P       ,
     &                     DA1M       ,DA2C       ,DA2P       ,
     &                     DA2M       ,SFNL       ,DSNL       ,
     &                     MEMNL4     ,WWINT      ,WWAWG      ,
     &                     WWSWG      ,CGO        ,USTAR      ,
     &                     ZELEN      ,SPCDIR     ,ANYWND     ,
     &                     DISSC0     ,DISSC1     ,GENC0      ,           40.85
     &                     GENC1      ,REDC0      ,REDC1      ,           40.85
     &                     XIS        ,FRCOEF     ,IT         ,           30.00
     &                     NPLA2      ,                                   40.55
     &                     URSELL     ,ANYBIN     ,REFLSO                 40.41 40.03
     &                                                        )           30.21
!
!****************************************************************

      USE OCPCOMM1                                                        40.41
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM3                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SWCOMM1                                                         40.41
      USE SWCOMM2                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_SNL4                                                          40.17

      IMPLICIT NONE                                                       40.17
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
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     32.06: Roeland Ris
!     40.02: IJsbrand Haagsma
!     40.03: Nico Booij
!     40.12: IJsbrand Haagsma
!     40.17: IJsbrand Haagsma
!     40.22: John Cazes and Tim Campbell
!     40.23: Marcel Zijlema
!     40.41: Andre van der Westhuysen
!     40.41: Marcel Zijlema
!     40.55: Marcel Zijlema
!     40.61: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     20.72, Jan. 96: Common introduced
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description several variables
!     32.06, June 99: Updated argument list of WNDPAR
!     30.81, Sep. 99: Updated argument list of SSURF
!     40.03, Apr. 00: array Ursell added in argument list
!     40.02, Sep. 00: Replaced SWCAP1-5 by SWCAP
!     40.02, Oct. 00: References to CDRAGP and TAUWP removed
!     40.12, Nov. 00: Added WCAP to dissipation output (bug fix 40.11 A)
!     40.22, Sep. 01: Removed WAREA array.                                40.22
!     40.22, Sep. 01: Changed array definitions to use the parameter      40.22
!                     MICMAX instead of ICMAX.                            40.22
!     40.17, Dec. 01: Implemented Multiple DIA                            40.17
!     40.23, Aug. 02: Print of CPU times added
!     40.23, Aug. 02: Parameter list of SWSNL2 and FILNL3 changed
!     40.41, May  04: Implementation of XNL (WRT) interface
!     40.41, Aug. 04: contribution due to reflection added to right-hand
!                     side of the system of equations
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.55, Dec. 05: introducing vegetation model
!     40.61, Sep. 06: introduction of all separate dissipation coefficients
!                     for output purposes
!     40.85, Aug. 08: add generation and redistribition for output purposes
!
!  2. Purpose
!
!     to compute the source terms, i.e., bottom friction,
!     wave breaking, wind input, white capping and non linear
!     wave wave interactions
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! i   ECOS  : =SPCDIR(*,2); cosine of spectral directions
! i   ESIN  : =SPCDIR(*,3); sine of spectral directions
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    ECOS(MDC)
      REAL    ESIN(MDC)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
!
!     INTEGERS:
!     --------------------------------------------------------------
!     IS       Counter of relative frequency band
!     IBOT     Indicator for bottom friction
!     ICUR     Indicator for current
!     ISURF    Indicator for wave breaking
!     ITRIAD   Indicator for nonlinear triad interactions
!     IQUAD    Indicator for nonlinear quadruplet interactions
!     IWCAP    Indicator for wave capping
!     IWIND    Indicator for which wind generation model is used
!     ICMAX    Maximum array size for the points of the molecul
!     MSC      Maximum counter of relative frequency in
!              computational model
!     MDC      Maximum counter of directional distribution in
!              computational model (2PI / DDIR + 1)
!     MTC      Maximum counter of the time, i.e.:
!              (total time in proto type) / (time step)
!     MBOT     Maximum array size for PBOT
!     MSURF    Maximum array size for PSSURF
!     MTRIAD   Maximum array size for PTRIAD
!     MWCAP    Maximum array size for PWCAP
!     MWIND    Maximum array size for PWIND
!     ISSTOP   Max frequency that is propagated within a sweep
!
!     REALS:
!     --------
!
!     DS          Width of frequency band (is not constant because
!                 of the logharitmic distribution of the frequency
!                 direction of the sweep (+1. OR -1. ) no input
!     GRAV        Gravitational acceleration
!     ABRBOT      Near bottom excursion amplitude
!     EMAX        Maximum energy according to the depth and the
!                 breaker parameter
!     ETOT        Total energy density per gridpoint
!     ETOTW       Total energy of the wind sea spectrum
!     GRAV        Gravitational acceleration
!     HM          Maximum wave height
!     KMESPC      Mean average wavenumber over full spectrum
!     SMESPC      Mean average frequency over full spectrum
!     QBLOC       Fraction of breaking waves
!
!     one and more dimensional arrays:
!     ---------------------------------
!
!     AC2       4D    (Nonstationary case) action density as function
!                     of D,S,X,Y at time T+DT
!     DEP2      2D    (Nonstationary case) depth as function of X and Y
!                     at time T+DIT
!     ECOS      1D    Represent the values of cos(d) of each spectral
!                     direction
!     ESIN      1D    Represent the values of sin(d) of each spectral
!                     direction
!     ALIMW     1D    Maximum energy by wind growth.
!     IMATDA    2D    coefficients of main diagonal
!     IMATRA    2D    right-hand side
!     KWAVE     2D    wavenumber as function of the relative frequency S
!                     and position IC(ix,iy)
!     PBOT      1D    Coefficient for the bottom friction models
!     PSURF     1D    Coefficient for the wave breaking model
!     PTRIAD    1D    Coefficient for the triad interaction model
!     PWCAP     1D    Coefficient for the white capping model
!     PWIND     1D    Coefficient for the wind growth model
!     UBOT      2D    Absolute orbital velocity in a gridpoint (IX,IY)
!     UX2       2D    (Nonstationary case) X-component of current velocity
!                     in (X,Y) at time T+DIT
!     UY2       2D    (Nonstationary case) Y-component of current velocity
!                     in (X,Y) at time T+DIT
!     USTAR     2D    Friction velocity at previous iteration for
!                     Janssen (1989,1990) wind input formulation
!     ZELEN     2D    Roughness length at previous iteration for
!                     Janssen (1989,1990) wind input formulation
!
!     Coefficients for the arrays:
!     -----------------------------
!                         default
!                         value:
!
!     PBOT(1)   = CFC      0.005    (Putnam and Collins equation)
!     PBOT(2)   = CFW      0.01     (Putnam and Collins equation)
!     PBOT(3)   = GAMJNS   0.067    (Jonswap equation)
!     PBOT(4)   = MF      -0.08     (Madsen et al. equation)
!     PBOT(5)   = KN       0.05     (Madsen et al. bottom roughness)
!
!     PSURF(1)  = ALFA     1.0      (Battjes & Janssen, 1978)
!     PSURF(2)  = GAMMA    0.73     (breaking criterium)
!
!     PWCAP(1)  = ALFAWC   2.36e-5  (Emperical coefficient)
!     PWCAP(2)  = ALFAPM   3.02E-3  (Alpha of Pierson Moskowitz frequency)
!     PWCAP(3)  = CFJANS   4.5
!     PWCAP(4)  = DELTA    0.5
!     PWCAP(5)  = CFLHIG   1.
!     PWCAP(6)  = GAMBTJ   0.88     (Steepness limited wave breaking )
!
!     PWIND(1)  = CF10     188.0    (second generation wind growth model)
!     PWIND(2)  = CF20     0.59     (second generation wind growth model)
!     PWIND(3)  = CF30     0.12     (second generation wind growth model)
!     PWIND(4)  = CF40     250.0    (second generation wind growth model)
!     PWIND(5)  = CF50     0.0023   (second generation wind growth model)
!     PWIND(6)  = CF60    -0.2233   (second generation wind growth model)
!     PWIND(7)  = CF70     0.       (second generation wind growth model)
!     PWIND(8)  = CF80    -0.56     (second generation wind growth model)
!     PWIND(9)  = RHOAW    0.00125  (density air / density water)
!     PWIND(10) = EDMLPM   0.0036   (limit energy Pierson Moskowitz)
!     PWIND(11) = CDRAG    0.0012   (drag coefficient)
!     PWIND(12) = UMIN     1.0      (minimum wind velocity)
!     PWIND(13) = PMLM     0.13     (  )
!
!     arrays for Janssen (`89)
!     -----------
!     PWIND(14) 1D    alfa (which is tuned at 0.01)
!     PWIND(15) 1D    Kappa ( 0.41)
!     PWIND(16) 1D    Rho air (1.28)
!     PWIND(17) 1D    Rho water (1025)
!
!  6. Local variables
!
!     IT    : Number of the time-step                                     40.17
!
!     IQERR : Error indicator for SWINTFXNL interface                     40.41
      INTEGER IQERR                                                       40.41
!
      INTEGER           :: ID, IDIA, IDC, IERR, IENT, IS, ISC, IT         40.17
      INTEGER           :: N2, LMAX                                       40.17

      REAL              :: DQ, DQ2, DT2                                   40.17
!
!  7. Common blocks used
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   If SBOT is on (IBOT > 0 ) then,
!     Call SBOT  to compute the source term due to bottom friction
!      according to Hasselmann et al. (1974), Putnam and Jonsson (1949)
!      or Madsen et al. (1991)
!   ------------------------------------------------------------
!   If SVEG is on (IVEG > 0 ) then,
!     Call SVEG  to compute the source term due to vegetation
!      dissipation according to Dalrymple (1984)
!   ------------------------------------------------------------
!   If SSURF is on (ISURF > 0 ) then,
!     Call SSURF to compute the source term due to wave breaking
!   ------------------------------------------------------------
!   IF IWIND =1 OR IWIND =2 THEN
!     Call WNDPAR (first or second generation mode of source terms
!                  using the DOLPHIN-B formulations)
!
!   else if IWIND = 3 then
!     input source term according to Snyder (1981)
!     Call SWIND3
!   else if IWIND = 4 then
!     input source term according to Janssen (1989,1991)
!     Call SWIND4
!   else if IWIND = 5 then
!     input source term according to Yan (1989) [reduces to Snyder form
!     for low frequencies and to Plant's (1982) form for high freq.
!     Call SWIND5
!   ------------------------------------------------------------
!   If IWCAP > 1 then
!     Call SWCAP to compute the source term for white capping             40.02
!   ---------------------------------------------------------------------
!   If ITRIAD > 0 then                                                    30.80
!      Call SWLTA to compute the nonlinear 3 wave-wave interactions
!      based on the LTA technique (Eldeberky, 1996)
!   ---------------------------------------------------------------------
!   If Ursell < Urmax
!   Then If IQUAD = 1
!        Then Call SWSNL1 to compute the nonlinear 4-wave interactions
!             semi-implicit per sweep direction
!        Else if IQUAD = 2
!        Then Call SWSNL2 to compute the nonlinear 4-wave interactions
!             fully explicit per sweep direction
!        Else if IQUAD = 3
!        Then Call SWSNL3 to compute the nonlinear 4-wave interactions
!             fully explicit per iteration
!             Call FILNL3 to get values for interactions from array
!             for full circle
!   ---------------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  ITER    ,IDWMIN  ,IDWMAX  ,SWPDIR  ,ISSTOP  ,
     &         IDDTOP  ,IDDLOW  ,IX      ,IY
!
      REAL     ABRBOT  ,ETOT    ,HM      ,QBLOC   ,ETOTW   ,
     &         FPM     ,WIND10  ,THETAW  ,SMESPC  ,KMESPC  ,
     &         SNLC1   ,FACHFR  ,DAL1    ,DAL2    ,DAL3    ,
     &         UFRIC   ,SMEBRK  ,HS      ,XIS
!
!
      REAL  :: AC2(MDC,MSC,MCGRD)                                         30.21
      REAL  :: DEP2(MCGRD)
      REAL  :: ALIMW(MDC,MSC)
      REAL  :: IMATDA(MDC,MSC)
      REAL  :: IMATRA(MDC,MSC)
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL  :: KWAVE(MSC,MICMAX)                                          40.22
      REAL  :: UBOT(MCGRD)
      REAL  :: UX2(MCGRD)
      REAL  :: UY2(MCGRD)
      REAL  :: UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )
      REAL  :: MEMNL4(MDC,MSC,MCGRD)
      REAL  :: PLWNDA(MDC,MSC,NPTST)
      REAL  :: PLWNDB(MDC,MSC,NPTST)
      REAL  :: PLWCAP(MDC,MSC,NPTST)
      REAL  :: PLBTFR(MDC,MSC,NPTST)
      REAL  :: PLWBRK(MDC,MSC,NPTST)
      REAL  :: PLNL4S(MDC,MSC,NPTST)
      REAL  :: PLNL4D(MDC,MSC,NPTST)
      REAL  :: PLVEGT(MDC,MSC,NPTST)                                      40.55
      REAL  :: PLTRI(MDC,MSC,NPTST)
      REAL  :: WWAWG(*)
      REAL  :: WWSWG(*)
      REAL  :: CGO(MSC,MICMAX)                                            40.22
      REAL  :: USTAR(MCGRD)
      REAL  :: ZELEN(MCGRD)
      REAL  :: DISSC0(1:MDC,1:MSC,1:MDISP)                                40.67
      REAL  :: DISSC1(1:MDC,1:MSC,1:MDISP)                                40.67
      REAL  :: GENC0 (1:MDC,1:MSC,1:MGENR)                                40.85
      REAL  :: GENC1 (1:MDC,1:MSC,1:MGENR)                                40.85
      REAL  :: REDC0 (1:MDC,1:MSC,1:MREDS)                                40.85
      REAL  :: REDC1 (1:MDC,1:MSC,1:MREDS)                                40.85
      REAL  :: URSELL(MCGRD)                                              40.03
      REAL  :: FRCOEF(MCGRD)                                              20.68
      REAL  :: NPLA2(MCGRD)                                               40.55
      REAL  :: REFLSO(MDC,MSC)                                            40.41
!

      INTEGER  IDCMIN(MSC)    ,
     &         IDCMAX(MSC)    ,
     &         WWINT(*)
!
      LOGICAL  GROWW(MDC,MSC) ,
     &         ANYBIN(MDC,MSC),                                           40.41
     &         ANYWND(MDC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SOURCE')
!
!     *** set relevant matrix elements to 0 ***
!
      IF (OPTG.NE.5) THEN                                                 40.80
         IMATRA = 0.
         IMATDA = 0.
      ENDIF
!
!     *** set all dissipation coeff at 0 ***
!
      DISSC0(1:MDC,1:MSC,1:MDISP) = 0.                                    40.67
      DISSC1(1:MDC,1:MSC,1:MDISP) = 0.                                    40.67
!
!     *** set all generation coeff at 0 ***
!
      GENC0(1:MDC,1:MSC,1:MGENR) = 0.                                     40.85
      GENC1(1:MDC,1:MSC,1:MGENR) = 0.                                     40.85
!
!     *** set all redistribution coeff at 0 ***
!
      REDC0(1:MDC,1:MSC,1:MREDS) = 0.                                     40.85
      REDC1(1:MDC,1:MSC,1:MREDS) = 0.                                     40.85
!
!TIMG      CALL SWTSTA(130)                                                    40.23
      IF (IBOT .GE. 1) THEN
!
!       *** wave-bottom interactions ***
!
        CALL SBOT (ABRBOT   ,DEP2     ,ECOS     ,ESIN     ,AC2      ,     41.04
     &             IMATDA   ,KWAVE    ,SPCSIG   ,UBOT     ,UX2      ,     30.72
     &             UY2      ,IDCMIN   ,IDCMAX   ,
     &             PLBTFR   ,ISSTOP   ,DISSC1   ,VARFR    ,FRCOEF   )     40.67
      END IF
!TIMG      CALL SWTSTO(130)                                                    40.23
!
!TIMG      CALL SWTSTA(139)                                                    40.55
      IF ( IVEG.GE.1 ) THEN                                               40.55
!
!     *** wave-vegetation interactions ***                                40.55
!
!        *** energy dissipation according to Dalrymple (1984)             40.55
!
         CALL SVEG (DEP2   ,IMATDA   ,ETOT   ,SMEBRK    ,                 40.58 40.55
     &              KMESPC ,PLVEGT   ,                                    40.58 40.55
     &              IDCMIN ,IDCMAX   ,ISSTOP ,DISSC1    ,                 40.55
     &              NPLA2  )                                              40.67 40.61 40.55
      END IF                                                              40.55
!TIMG      CALL SWTSTO(139)                                                    40.55
!
!TIMG      CALL SWTSTA(131)                                                    40.23
      IF (ISURF .GE. 1) THEN
!
!         *** calculate surf breaking source term (5 formulations) ***    41.03
!
          CALL SSURF (ETOT    ,HM      ,QBLOC   ,SMEBRK  ,
     &                KMESPC  ,SPCSIG  ,AC2     ,IMATRA  ,                30.81
     &                IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,
     &                ISSTOP  ,DISSC0  ,DISSC1  )                         40.67 40.61 30.21
!
      END IF
!TIMG      CALL SWTSTO(131)                                                    40.23
!
!TIMG      CALL SWTSTA(132)                                                    40.23
      IF ( IWIND .GE. 3
     &                 ) THEN
!
!       *** linear wind input according to Cavaleri and Malanotte ***
!       *** Rizolli (1981) for a third generation mode of SWAN    ***
!
        IF (PWIND(31) .GT. 1.E-20)                                        7/MAR
     &  CALL SWIND0 (IDCMIN  ,IDCMAX  ,ISSTOP  ,
     &               SPCSIG  ,THETAW  ,ANYWND  ,
     &               UFRIC   ,FPM     ,PLWNDA  ,
     &               IMATRA  ,SPCDIR  ,GENC0   )
      ENDIF
!
      IF ( IWIND .EQ. 1 .OR. IWIND .EQ. 2 ) THEN                          970220
!
         CALL WNDPAR (ISSTOP,IDWMIN,IDWMAX,IDCMIN,IDCMAX,                 32.06
     &                   DEP2  ,WIND10,GENC0,GENC1,                       40.85 32.06
     &                   THETAW,AC2   ,KWAVE ,IMATRA,IMATDA,              32.06
     &                   SPCSIG,CGO   ,ALIMW ,GROWW ,ETOTW ,              32.06
     &                   PLWNDA,PLWNDB,SPCDIR,ITER            )           32.06
!
!
      ELSE IF ( IWIND .EQ. 3 ) THEN
!
!       *** Wind input according to Snyder et al (1981) ***
!
        CALL SWIND3 (SPCSIG  ,THETAW  ,
     &               KWAVE   ,IMATRA  ,GENC0   ,
     &               IDCMIN  ,IDCMAX  ,AC2     ,UFRIC   ,
     &               FPM     ,PLWNDA  ,ISSTOP  ,SPCDIR  ,ANYWND  )
!
      ELSE IF ( IWIND .EQ. 4 ) THEN
!
!       *** Wind input according to Janssen (1989,1991) ***
!
        CALL SWIND4  (IDWMIN  ,IDWMAX  ,
     &                SPCSIG  ,WIND10  ,THETAW  ,XIS     ,
     &                DDIR    ,KWAVE   ,IMATRA  ,GENC0   ,
     &                IDCMIN  ,IDCMAX  ,AC2     ,UFRIC   ,
     &                PLWNDA  ,ISSTOP  ,ITER    ,USTAR   ,ZELEN   ,
     &                SPCDIR  ,ANYWND  ,IT      )
!
      ELSE IF ( IWIND .EQ. 5 ) THEN
!
!       *** Wind input according to Yan (1989) ***
!
        CALL SWIND5 (SPCSIG  ,THETAW  ,ISSTOP  ,
     &               UFRIC   ,KWAVE   ,IMATRA  ,IDCMIN  ,
     &               IDCMAX  ,AC2     ,ANYWND  ,PLWNDA  ,
     &               SPCDIR  ,GENC0                     )

      END IF
!TIMG      CALL SWTSTO(132)                                                    40.23
!
!     Calculate whitecapping source term (six formulations)               40.31 40.02
!
!TIMG      CALL SWTSTA(133)                                                    40.23
      IF (IWCAP.GE.1) CALL SWCAP (SPCDIR  ,SPCSIG  ,KWAVE   ,AC2     ,    40.02
     &                            IDCMIN  ,IDCMAX  ,ISSTOP  ,             40.02
     &                            ETOT    ,IMATDA  ,IMATRA  ,PLWCAP  ,    40.02
     &                            CGO     ,UFRIC   ,                      40.53
     &                            DEP2    ,DISSC1  ,DISSC0  )             40.67 40.61 40.12
!TIMG      CALL SWTSTO(133)                                                    40.23
!
!     compute nonlinear interactions, starting with triads                 NB!
!
!TIMG      CALL SWTSTA(134)                                                    40.23
      IF (ITRIAD .GT. 0) THEN
!
!       *** compute the 3 wave-wave interactions if in each ***            NB!
!       *** geographical gridpoint a continuous spectrum    ***
!       *** is present, i.e., after first iteration         ***
!
        IF ( ICUR .EQ. 0 .OR. ITER .GT. 1 ) THEN
!
            CALL SWLTA ( AC2   , DEP2  , CGO   , SPCSIG,                  40.56
     &                   KWAVE , IMATRA, IMATDA, REDC0 , REDC1 ,          40.85
     &                   IDDLOW, IDDTOP, ISSTOP, IDCMIN, IDCMAX,
     &                   HS    , SMEBRK, PLTRI , URSELL )
!
        ENDIF
!
      ENDIF
!TIMG      CALL SWTSTO(134)                                                    40.23
!
!     --- compute quadruplet interactions if Ursell number < Urmax        40.41
!         (usually, Urmax = 0.1, but here Urmax = 10)                     40.41
!
!TIMG      CALL SWTSTA(135)                                                    40.23
      IF (URSELL(KCGRD(1)).LT.PTRIAD(3)) THEN                             40.03
!
!       *** compute the counters for the nonlinear four ***
!       *** wave-wave interactions in spectral space    ***
!       *** and high frequency factor                   ***
!
        IF ( IQUAD .GE. 1 ) THEN
           CALL RANGE4 (WWINT ,IDDLOW,IDDTOP )                            40.00
           FACHFR = 1. / XIS ** PWTAIL(1)                                 20.72
        ENDIF

!
        IF (IQUAD .EQ. 1) THEN
!
!       *** semi-implicit calculation for al the bins that fall ***
!       *** within a sweep. No additional array is required     ***
!
          CALL SWSNL1 (                  WWINT   ,WWAWG   ,WWSWG   ,      34.00
     &                 IDCMIN  ,IDCMAX  ,UE      ,SA1     ,               40.17
     &                 SA2     ,DA1C    ,DA1P    ,DA1M    ,DA2C    ,
     &                 DA2P    ,DA2M    ,SPCSIG  ,SNLC1   ,KMESPC  ,      30.72
     &                 FACHFR  ,ISSTOP  ,DAL1    ,DAL2    ,DAL3    ,
     &                 SFNL    ,DSNL    ,DEP2    ,AC2     ,IMATDA  ,
     &                 IMATRA  ,PLNL4S  ,PLNL4D                    ,      34.00
     &                 IDDLOW  ,IDDTOP  ,REDC0   ,REDC1   )               40.85 34.00
!
        ELSE IF ( IQUAD .EQ. 2) THEN
!
!         *** fully explicit calculation for al the bins that fall ***
!         *** within a sweep. No additional array is required      ***
!
          CALL SWSNL2 (                  IDDLOW  ,IDDTOP  ,WWINT   ,      34.00
     &                 WWAWG   ,UE      ,SA1     ,ISSTOP  ,               40.17
     &                 SA2     ,SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,      30.72
     &                 DAL3    ,SFNL    ,DEP2    ,AC2     ,KMESPC  ,
     &                          REDC0   ,REDC1   ,IMATDA  ,IMATRA  ,      40.85 40.23 34.00
     &                 FACHFR  ,PLNL4S  ,         IDCMIN  ,IDCMAX  )      34.00
!
        ELSE IF ( IQUAD .EQ. 3) THEN
!
!         *** fully explicit calculation of the 4 wave-wave inter-  ***
!         *** actions for the full circle (1 -> MDC). An additional ***
!         *** array is required in which the values are stored prior***
!         *** to every iteration                                    ***
!
          IF ( ITER .EQ. 1 ) THEN
!
!           *** calculate the interactions every sweep in each grid ***
!           *** point for the first iteration to ensure stable      ***
!           *** behaviour of the model                              ***
!
            CALL SWSNL3 (                  WWINT   ,WWAWG   ,             40.17
     &                   UE      ,SA1     ,SA2     ,SPCSIG  ,SNLC1   ,    40.17 30.72
     &                   DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    ,    40.17
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR           )    40.17
!
          ELSE IF ( ITER .GT. 1 .AND. ( SWPDIR .EQ. 1 .OR.
     &        ( SWPDIR .EQ. 2 .AND. IX .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 3 .AND. IY .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 4 .AND. (IX.EQ.MXC .AND. IY.EQ.1)) )) THEN
!
            CALL SWSNL3 (                  WWINT   ,WWAWG   ,             40.17
     &                   UE      ,SA1     ,SA2     ,SPCSIG  ,SNLC1   ,    40.17 30.72
     &                   DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    ,    40.17
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR           )    40.17
!
          ENDIF
!
!         *** Get source term value of additional array for the bin   ***
!         *** that fall within a sweep and store in right hand vector ***
!
          CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,AC2     ,      40.23
     &                 MEMNL4  ,PLNL4S  ,ISSTOP  ,REDC0   ,REDC1   )      40.85 40.41

        ELSE IF (IQUAD .EQ. 4) THEN                                       40.17

!         Multiple DIA according to Hashimoto (1999)                      40.17

          IF ((ITER.EQ.1).OR.( ITER .GT. 1 .AND. ( SWPDIR .EQ. 1 .OR.
     &        ( SWPDIR .EQ. 2 .AND. IX .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 3 .AND. IY .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 4 .AND. (IX.EQ.MXC .AND. IY.EQ.1))))) THEN
            DO IDIA=1,MDIA                                                40.17
              PQUAD(1) = LAMBDA(IDIA)                                     40.17
              CALL FAC4WW (XIS   ,SNLC1 ,                                 40.41 40.17
     &                     DAL1  ,DAL2  ,DAL3  ,SPCSIG,                   40.17
     &                     WWINT ,WWAWG ,WWSWG                )           40.17
              FACHFR = 1. / XIS ** PWTAIL(1)                              40.17
              CALL RANGE4 (WWINT ,IDDLOW,IDDTOP )                         40.17
              CALL SWSNL4 (WWINT   ,WWAWG   ,                             40.17
     &                     SPCSIG  ,SNLC1   ,                             40.17
     &                     DAL1    ,DAL2    ,DAL3    ,DEP2    ,           40.17
     &                     AC2     ,KMESPC  ,MEMNL4  ,FACHFR  ,           40.17
     &                     IDIA    ,ITER )                                40.17
            END DO
          ENDIF

!         Fill the matrix per sweep even though the quadruplets are calculated
!         only once per iteration

          CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,AC2     ,      40.23
     &                 MEMNL4  ,PLNL4S  ,ISSTOP  ,REDC0   ,REDC1   )      40.85 40.41
!
        ELSE IF ( IQUAD .EQ. 8) THEN                                      40.41
!
!         --- fully explicit calculation of the 4 wave-wave
!             interactions for the full circle. The interactions
!             in neighbouring bins are interpolated in piecewise
!             constant manner. An additional array is required in
!             which the values are stored prior to every iteration
!
          IF ( ITER .EQ. 1 ) THEN
!
!           *** calculate the interactions every sweep in each grid ***
!           *** point for the first iteration to ensure stable      ***
!           *** behaviour of the model                              ***
!
            CALL SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,
     &                   SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,
     &                   DEP2    ,AC2     ,KMESPC  ,MEMNL4  ,FACHFR  )
!
          ELSE IF ( ITER .GT. 1 .AND. ( SWPDIR .EQ. 1 .OR.
     &        ( SWPDIR .EQ. 2 .AND. IX .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 3 .AND. IY .EQ. 1) .OR.
     &        ( SWPDIR .EQ. 4 .AND. (IX.EQ.MXC .AND. IY.EQ.1)) )) THEN
!
            CALL SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,
     &                   SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,
     &                   DEP2    ,AC2     ,KMESPC  ,MEMNL4  ,FACHFR  )
!
          ENDIF
!
!         *** get source term value of additional array for the bin   ***
!         *** that fall within a sweep and store in right hand vector ***
!
          CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,AC2     ,
     &                 MEMNL4  ,PLNL4S  ,ISSTOP  ,REDC0   ,REDC1   )      40.85 40.41
!
        ELSEIF ((IQUAD.EQ.51).OR.(IQUAD.EQ.52).OR.(IQUAD.EQ.53)) THEN     40.41
!
!         Calculate the quadruplets using the XNL interface of            40.41
!         G. van Vledder                                                  40.41
!
!         Avoid calculation of the quadruplets in more than one sweep     40.41
!
          IF ((ITER .GE. 1) .AND.                                         40.41
     &         ( (SWPDIR.EQ.1)                                 .OR.       40.41
     &          ((SWPDIR.EQ.2).AND.(IX.EQ.1)                  ).OR.       40.41
     &          ((SWPDIR.EQ.3).AND.(IY.EQ.1)  .AND.(.NOT.ONED)).OR.       40.41
     &          ((SWPDIR.EQ.4).AND.(IX.EQ.MXC).AND.(IY.EQ.1)              40.41
     &                                        .AND.(.NOT.ONED)) )         40.41
     &       ) THEN                                                       40.41
!
            IF (ITEST.GE.30) WRITE(PRINTF,'(A,4I6,F12.2)')                40.41
     &         'SOURCE XNL: iter, swpdir, kcgrd(1), iquad, depth:',       40.41
     &         ITER, SWPDIR, KCGRD(1), IQUAD, DEP2(KCGRD(1))              40.41

            CALL SWINTFXNL(AC2,SPCSIG,SPCDIR,MDC,MSC,MCGRD,               40.41
     &                     DEP2,IQUAD,MEMNL4,KCGRD,ICMAX,IQERR)           40.41
!
          ENDIF                                                           40.41
!
          IF (ITEST.GE.30) THEN                                           40.41
             WRITE (PRTEST,*) '+SOURCE: IX, IY, SWPDIR: ',                40.41
     &                         IX, IY, SWPDIR                             40.41
          ENDIF                                                           40.41
          IF (TESTFL.AND.ITEST.GE.100) THEN                               40.41
             DO IS=1, MSC                                                 40.41
                DO ID = 1, MDC                                            40.41
                  WRITE(PRINTF,9300) IS,ID,MEMNL4(ID,IS,KCGRD(1))         40.41
 9300             FORMAT(' SOURCE: IS ID MEMNL(): ',2I6,E12.4)            40.41
                ENDDO                                                     40.41
             ENDDO                                                        40.41
          ENDIF                                                           40.41
!
          CALL FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,AC2     ,      40.41
     &                 MEMNL4  ,PLNL4S  ,ISSTOP  ,REDC0   ,REDC1   )      40.85 40.41
!
          IF (ITEST.GE.30) THEN                                           40.41
            WRITE (PRTEST,*) '+SOURCE: ITER, IQUAD, SWPDIR, IQERR: ',     40.41
     &                       ITER, IQUAD, SWPDIR, IQERR                   40.41
          ENDIF                                                           40.41
!
        ENDIF
      ENDIF
!TIMG      CALL SWTSTO(135)                                                    40.23
!
!     --- add contribution due to reflection of obstacles
!
!TIMG      CALL SWTSTA(136)                                                    40.41
      IF (NUMOBS.NE.0) THEN                                               40.41
         DO IS = 1, MSC                                                   40.41
            DO ID = 1, MDC                                                40.41
               IF (ANYBIN(ID,IS))                                         40.41
     &            IMATRA(ID,IS) = IMATRA(ID,IS) + REFLSO(ID,IS)           40.41
            END DO                                                        40.41
         END DO                                                           40.41
      END IF                                                              40.41
!TIMG      CALL SWTSTO(136)                                                    40.41
!
!     End of the subroutine SOURCE
      RETURN
      END
!
!************************************************************************
!                                                                       *
      SUBROUTINE PHILIM(AC2,AC2OLD,CGO,KWAVE,SPCSIG,ANYBIN,ISLMIN,NFLIM,
     &                  QB_LOC)
!                                                                       *
!************************************************************************
!
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
!     30.82: IJsbrand Haagsma
!     40.16: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.82, Feb. 99: New subroutine
!     40.16, Dec. 01: Implemented limiter switch
!     40.23, Aug. 02: Store number of frequency use of limiter
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Limits the change in action density between two iterations to a
!     certain percentage of the (directionally independent) Phillips
!     equilibrium level
!
!  3. Method
!
!     The maximum change of energy density per bin is related to
!     the (directionally independent) Phillips equilibrium level.
!     This change is estimated as:
!
!     |D E(s)| = factor * alpha_PM * g^2 / (s^5)
!
!     in which the Phillips' constant for a Pierson-Moskowitz
!     spectrum (alpha_PM) is taken to be 0.0081. Note that this
!     is a measure for a 1D spectrum. In SWAN, factor = 0.1
!     (stored in PNUMS(20)).
!
!     In terms of action density, we have
!
!     |D N(s)| = factor * alpha_PM * g^2 / (s^6)
!
!     Expressing in wave number k, this becomes with a deep
!     water approach of s^2 = gk:
!
!     |D N(s)| = factor * alpha_PM / (s^2 k^2)
!
!     Furthermore, with s = 2*k*c_g (deep water), we finally
!     have (Ris, 1997, p.36):
!
!                          alpha_PM
!     |D N(s,t)| = factor -----------
!                         2 s k^3 c_g
!
!     In cases where waves are breaking the dissipation of energy
!     is not limited. This is assumed to be the case when the
!     fraction of breaking waves Qb is more than 1.e-5.
!
!  4. Argument variables
!
      LOGICAL ANYBIN(MDC,MSC)
!
!     ISLMIN: Lowest sigma-index occured in applying limiter
!     NFLIM : Number of frequency use of limiter
!     QB_LOC: Local value of Qb (fraction of breaking waves)
!
      INTEGER ISLMIN(MCGRD), NFLIM(MCGRD)
      REAL    QB_LOC
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    AC2OLD(MDC,MSC)
      REAL    CGO(MSC,ICMAX)
      REAL    KWAVE(MSC,ICMAX)
      REAL    SPCSIG(MSC)
!
!  6. Local variables
!
!     ID    : Counter for directional (theta) space
!     IS    : Counter for frequency (sigma) space
!
      INTEGER ID,IS
!
!     DAC2MX: Maximum deviation of action density AC2 between iterations
!
      REAL    DAC2MX
!
! 13. Source text
!
      IF (MSC.GT.3) THEN
         DO IS=1,MSC
            DAC2MX=ABS((PNUMS(20)*0.0081)/
     &             (2.*SPCSIG(IS)*(KWAVE(IS,1)**3)*CGO(IS,1)))
            DO ID=1,MDC
               IF (ANYBIN(ID,IS) .AND.
     &             AC2(ID,IS,KCGRD(1)).GT.AC2OLD(ID,IS)+DAC2MX) THEN
                  AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)+DAC2MX
                  NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                  ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
               END IF
            END DO
            IF (QB_LOC.LT.PNUMS(28)) THEN                                 40.16
               DO ID=1,MDC
                  IF (ANYBIN(ID,IS) .AND.
     &                AC2(ID,IS,KCGRD(1)).LT.AC2OLD(ID,IS)-DAC2MX) THEN
                     AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)-DAC2MX
                     NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                     ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
                  END IF
               END DO
            END IF
         END DO
      END IF
      RETURN
      END
!************************************************************************
!                                                                       *
      SUBROUTINE HJLIM(AC2,AC2OLD,CGO,KWAVE,SPCSIG,ANYBIN,ISLMIN,NFLIM,
     &                 QB_LOC,USTAR)
!                                                                       *
!************************************************************************
!
      USE SWCOMM3
      USE TIMECOMM
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
!     40.61: Roop Lalbeharry
!
!  1. Updates
!
!     40.61, Nov. 06: New subroutine
!
!  2. Purpose
!
!     Limits the change in action density between two iterations to
!     the (directionally independent) Hersbach and Janssen (1999) limiter
!
!  3. Method
!
!     The maximum change of energy density per bin is related to
!     the (directionally independent) Hersbach and Janssen (1999) limiter
!     This change is estimated in terms of frequency and energy density as:
!
!     |D E(f,t)| = 3.0 * 1.0E-7 * g * u* * f_c * dt/ (f^4)
!
!     in which f_c is the model's cut off frequency and u* the
!     friction velocity and dt integration time step (sec)
!
!     In terms of action density and angular frequency and for deep water,
!     we have:
!
!     |D N(s,t)| = C_HJ * u* / (s^3 * k)
!
!     where C_HJ = (2.0 * PI)^2 * 3.0 * 1.0E-7 * s_c * dt and
!     u* = max(u*, g*s*_pm/s); s*_pm = 2.*PI*f*_pm; f*_pm = 5.6 * 1.0E-3
!
!     Compare with Ris's(1997, p.36) formulation for deep water:
!
!                          alpha_PM
!     |D N(s,t)| = factor ----------- ; factor = 0.1, alpha_PM = 0.0081
!                         2 s k^3 c_g
!
!     In cases where waves are breaking the dissipation of energy
!     is not limited. This is assumed to be the case when the
!     fraction of breaking waves Qb is more than 1.e-5.
!
!  4. Argument variables
!
      LOGICAL ANYBIN(MDC,MSC)
!
!     ISLMIN: Lowest sigma-index occured in applying limiter
!     NFLIM : Number of frequency use of limiter
!     QB_LOC: Local value of Qb (fraction of breaking waves)
!
      INTEGER ISLMIN(MCGRD), NFLIM(MCGRD)
      REAL    QB_LOC
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    AC2OLD(MDC,MSC)
      REAL    CGO(MSC,ICMAX)
      REAL    KWAVE(MSC,ICMAX)
      REAL    SPCSIG(MSC)
      REAL    USTAR(MCGRD)
!
!  6. Local variables
!
!     ID    : Counter for directional (theta) space
!     IS    : Counter for frequency (sigma) space
!
      INTEGER ID,IS
!
!     DAC2MX: Maximum deviation of action density AC2 between iterations
!
      REAL    DAC2MX, UFRIC_VEL, C_HJ, SPM_NOND
!
! 13. Source text
!
      C_HJ = PI2**2*3.0*1.0E-7*DT*SPCSIG(MSC)
      SPM_NOND = PI2 * 5.6 * 1.0E-3
      IF (MSC.GT.3) THEN
         DO IS=1,MSC
            UFRIC_VEL = MAX(USTAR(KCGRD(1)), GRAV*SPM_NOND/SPCSIG(IS))
            DAC2MX=ABS((C_HJ*UFRIC_VEL)/
     &             (SPCSIG(IS)**3*KWAVE(IS,1)))
            DO ID=1,MDC
               IF (ANYBIN(ID,IS) .AND.
     &             AC2(ID,IS,KCGRD(1)).GT.AC2OLD(ID,IS)+DAC2MX) THEN
                  AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)+DAC2MX
                  NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                  ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
               END IF
            END DO
            IF (QB_LOC.LT.PNUMS(28)) THEN
               DO ID=1,MDC
                  IF (ANYBIN(ID,IS) .AND.
     &                AC2(ID,IS,KCGRD(1)).LT.AC2OLD(ID,IS)-DAC2MX) THEN
                     AC2(ID,IS,KCGRD(1))=AC2OLD(ID,IS)-DAC2MX
                     NFLIM(KCGRD(1)) = NFLIM(KCGRD(1)) + 1
                     ISLMIN(KCGRD(1)) = MIN(IS,ISLMIN(KCGRD(1)))
                  END IF
               END DO
            END IF
         END DO
      END IF
      RETURN
      END
!****************************************************************
!
      SUBROUTINE RESCALE (AC2, ISSTOP, IDCMIN, IDCMAX, NRSCAL)
!
!****************************************************************

      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
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
!  0. Authors
!
!     40.00: Nico Booij
!     40.23: Marcel Zijlema
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New subroutine (software moved from subroutines
!                     SOLBAND, SOLMT1 and SOLMAT
!     40.23, Aug. 02: Store number of frequency use of rescaling
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Remove negative values from a computed action density spectrum
!
!  3. Method
!
!     Make negative action densities 0 at the expense of other action densities
!     for the frequency
!
!  4. Argument variables
!
!     AC2         action densities
!
      REAL        AC2(MDC,MSC,MCGRD)
!
!     ISSTOP      maximum frequency counter in this sweep
!
      INTEGER     ISSTOP
!
!     IDCMIN      Integer array containing minimum counter of directions
!     IDCMAX      Integer array containing maximum counter
!     NRSCAL      Number of frequency use of rescaling
!
      INTEGER     IDCMIN(MSC), IDCMAX(MSC)
      INTEGER     NRSCAL(MCGRD)
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
!        ---
!
!     9. STRUCTURE
!
!   -------------------------------------------------------------
!   For all frequencies do
!       Make ATOT equal to integral of action density over direction
!       Make ATOTP equal to integral of positive action density
!       Determine FACTOR
!       If negative values do occur
!       Then for all directions do
!            If action density is negative
!            Then make action density =0
!            Else multiply action density by FACTOR
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
!         local variables
!
!         IS         counter of frequency
!         ID         counter of direction
!         IDDUM      uncorrected counter of direction
!
      INTEGER  IS      ,ID      ,IDDUM,   IENT
!
!         ATOT       integral of action density for one frequency
!         ATOTP      integral of positive action density for one frequency
!         FACTOR
!
      REAL     ATOT    ,ATOTP   ,FACTOR
!
!         NEGVAL      if True, there are negative values in the spectrum
!
      LOGICAL  NEGVAL
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'RESCALE')
!
!     *** if negative action density occur rescale with a factor ***
!     *** only the sector computed is rescaled !!                ***
!
      DO 180 IS = 1 , ISSTOP
        ATOT   = 0.
        ATOTP  = 0.
        FACTOR = 0.
        NEGVAL = .FALSE.
        DO 160 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          ATOT = ATOT + AC2(ID,IS,KCGRD(1))
          IF ( AC2(ID,IS,KCGRD(1)) .LT. 0. ) THEN
            NRSCAL(KCGRD(1)) = NRSCAL(KCGRD(1)) + 1
            NEGVAL = .TRUE.
          ELSE
            ATOTP = ATOTP + AC2(ID,IS,KCGRD(1))
          END IF
 160    CONTINUE
        IF (NEGVAL) THEN
          IF ( ATOTP .LT. 1.E-15 ) ATOTP = 1.E-15
          FACTOR = ATOT / ATOTP
!
!         *** rescale ***
!
          DO 170 IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( AC2(ID,IS,KCGRD(1)) .LT. 0.) THEN
              AC2(ID,IS,KCGRD(1)) = 0.
            END IF
            IF ( FACTOR .GE. 0. ) THEN
              AC2(ID,IS,KCGRD(1)) = FACTOR * AC2(ID,IS,KCGRD(1))
            ENDIF
 170      CONTINUE
!
          IF ( ITEST .GE. 120 .AND. TESTFL )
     &    WRITE (PRINTF, 171) IXCGRD(1)+MXF-2, IYCGRD(1)+MYF-2, IS,
     &    FACTOR , ATOT, ATOTP
 171      FORMAT(' Rescale in Point, Isig, Factor, ATOT, ATOTP:',
     &    3I4, 3(1X,E11.4))
        ENDIF
 180  CONTINUE
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSIP ( AC2   , IMATDA, IMATRA, IMATLA, IMATUA,
     &                   IMAT5L, IMAT6U, AC2OLD, REPS  , MAXIT ,
     &                   IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP,          40.41
     &                   IDCMIN, IDCMAX )                                 40.41
!
!****************************************************************
!
      USE SWCOMM1                                                         40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31

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
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Oct. 02: New subroutine
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.41, Mar. 04: parameter ALFA set to 0.0, extra test output
!                     and some corrections if SECTOR=0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Solves penta-diagonal system of equations in
!     spectral space by means of Stone's SIP solver
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     IAMOUT      control parameter indicating the amount of
!                 output required
!                 0: no output
!                 1: only fatal errors will be printed
!                 2: gives output concerning the iteration process
!                 3: additional information about the iteration
!                    is printed
!     IDCMAX      maximum counter in directional space
!     IDCMIN      minimum counter in directional space
!     IDDLOW      minimum direction that is propagated within a sweep
!     IDDTOP      maximum direction that is propagated within a sweep
!     IMAT5L      coefficients of lower diagonal in sigma-space
!     IMAT6U      coefficients of upper diagonal in sigma-space
!     IMATDA      coefficients of main diagonal
!     IMATLA      coefficients of lower diagonal in theta-space
!     IMATUA      coefficients of upper diagonal in theta-space
!     IMATRA      right-hand side
!     INOCNV      integer indicating number of grid points in which
!                 solver does not converged
!     ISSTOP      maximum frequency counter in a sweep
!     MAXIT       the maximum number of iterations to be performed in
!                 the linear solver
!     REPS        accuracy with respect to the right-hand side used
!                 in the following termination criterion:
!
!                 ||b-Ax || < reps*||b||
!                       k
!
      INTEGER IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP, MAXIT
      INTEGER IDCMIN(MSC), IDCMAX(MSC)
      REAL    REPS
      REAL    AC2(MDC,MSC,MCGRD),
     &        IMATDA(MDC,MSC), IMATRA(MDC,MSC),
     &        IMAT5L(MDC,MSC), IMAT6U(MDC,MSC),
     &        IMATLA(MDC,MSC), IMATUA(MDC,MSC),
     &        AC2OLD(MDC,MSC)
!
!  5. Parameter variables
!
!     ALFA        relaxation parameter used in the SIP solver
!     SMALL :     a small number
!
      REAL    ALFA, SMALL
      PARAMETER (ALFA=0.0,SMALL=1.E-15)                                   40.41
!
!  6. Local variables
!
!     BNORM :     2-norm of right-hand side vector
!     CMAT5L:     coefficients of lower diagonal in sigma-space
!                 obtained by an incomplete lower-upper factorization
!     CMAT6U:     coefficients of upper diagonal in sigma-space
!                 obtained by an incomplete lower-upper factorization
!     CMATDA:     coefficients of main diagonal obtained by an
!                 incomplete lower-upper factorization
!     CMATLA:     coefficients of lower diagonal in theta-space
!                 obtained by an incomplete lower-upper factorization
!     CMATUA:     coefficients of upper diagonal in theta-space
!                 obtained by an incomplete lower-upper factorization
!     EPSLIN:     required accuracy in the linear solver
!     ICONV :     indicator for convergence (1=yes, 0=no)
!     ID    :     loop counter
!     IDDL  :     minimum counter in theta-space of modulo MDC
!     IDDT  :     maximum counter in theta-space of modulo MDC
!     IDDUM :     loop counter
!     IDM   :     index of point ID-1
!     IDMAX :     local array of maximum counter in theta-space
!     IDMIN :     local array of minimum counter in theta-space
!     IDP   :     index of point ID+1
!     IENT  :     number of entries
!     IS    :     loop counter
!     ISM   :     index of point IS-1
!     ISP   :     index of point IS+1
!     IT    :     iteration count
!     LOPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!     P1    :     auxiliary factor
!     P2    :     auxiliary factor
!     P3    :     auxiliary factor
!     RES   :     the residual vector
!     RNORM :     2-norm of residual vector
!     RNRM0 :     2-norm of initial residual vector
!     UEPS  :     minimal accuracy based on machine precision
!     UPPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!
      INTEGER ICONV, ID, IDDL, IDDT, IDDUM, IDM, IDP, IENT,
     &        IS, ISM, ISP, IT
      INTEGER IDMIN(MSC), IDMAX(MSC)
      REAL    BNORM, EPSLIN, P1, P2, P3, RNORM, RNRM0, UEPS
      REAL    RES(MDC,MSC)   , CMATDA(MDC,MSC),
     &        CMAT5L(MDC,MSC), CMAT6U(MDC,MSC),
     &        CMATLA(MDC,MSC), CMATUA(MDC,MSC),
     &        LOPERI(MSC)    , UPPERI(MSC)
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
!     SWOMPU  (in SWANCOM1)
!
! 12. Structure
!
!     The system of equations is solved using an incomplete
!     factorization technique called Strongly Implicit Procedure
!     (SIP) as described in
!
!     H.L. Stone
!     Iterative solution of implicit approximations of
!     multidimensional partial differential equations
!     SIAM J. of Numer. Anal., vol. 5, 530-558, 1968
!
!     This method constructs an incomplete lower-upper factorization
!     that has the same sparsity as the original matrix. Hereby, a
!     parameter alfa is used, which should be 0.0 in case of SWAN         40.41
!     (when alfa > 0.95, the method may diverge).
!
!     Afterward, the resulting system is solved in an iterative manner
!     by forward and backward substitutions.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSIP')

!     --- initialize arrays

      RES    = 0.
      CMATDA = 0.
      CMATLA = 0.
      CMATUA = 0.
      CMAT5L = 0.
      CMAT6U = 0.
      LOPERI = 0.
      UPPERI = 0.

!     --- in case of periodicity in theta-space, store values
!         of matrix coefficients corresponding to left bottom and
!         right top

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).EQ.1 .AND. IDCMAX(IS).EQ.MDC ) THEN
           UPPERI(IS) = IMATLA(  1,IS)
           LOPERI(IS) = IMATUA(MDC,IS)
         END IF
      END DO

!     --- when no bins fall within the sweep, i.e. SECTOR = 0,
!         reset the bounds of sector as 1..MDC (routine SOLPRE
!         has clear the rows in the matrix that do not belong
!         to the sweep)

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).LE.IDCMAX(IS) ) THEN
            IDMIN(IS) = IDCMIN(IS)
            IDMAX(IS) = IDCMAX(IS)
         ELSE
            IDMIN(IS) = 1
            IDMAX(IS) = MDC
         END IF
      END DO

      IT    = 0
      ICONV = 0

!     --- construct L and U matrices (stored in CMAT[xx])

      BNORM = 0.

      IS    = 1
      IDDUM = IDMIN(IS)
      ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1

      CMAT5L(ID,IS) = IMAT5L(ID,IS)
      CMATLA(ID,IS) = IMATLA(ID,IS)
      CMATDA(ID,IS) = 1./(IMATDA(ID,IS)+SMALL)
      CMAT6U(ID,IS) = IMAT6U(ID,IS)*CMATDA(ID,IS)
      CMATUA(ID,IS) = IMATUA(ID,IS)*CMATDA(ID,IS)
      BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)

      DO IDDUM = IDMIN(IS)+1, IDMAX(IS)
         ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1

         P2 = ALFA*CMAT6U(IDM,IS)
         CMAT5L(ID,IS) = IMAT5L(ID,IS)
         CMATLA(ID,IS) = IMATLA(ID,IS)/(1.+P2)

         P2 = P2*CMATLA(ID,IS)
         P3 = IMATDA(ID,IS) + P2
     &       -CMATLA(ID,IS)*CMATUA(IDM,IS )
     &       +SMALL
         CMATDA(ID,IS) = 1./P3
         CMAT6U(ID,IS) = (IMAT6U(ID,IS)-P2)*CMATDA(ID,IS)
         CMATUA(ID,IS) =      IMATUA(ID,IS)*CMATDA(ID,IS)
         BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)
      END DO

      DO IS = 2, ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1

         P1 = ALFA*CMATUA(ID,ISM)
         CMAT5L(ID,IS) = IMAT5L(ID,IS)/(1.+P1)
         CMATLA(ID,IS) = IMATLA(ID,IS)
         P1 = P1*CMAT5L(ID,IS)
         P3 = IMATDA(ID,IS) + P1
     &       -CMAT5L(ID,IS)*CMAT6U(ID,ISM)
     &       +SMALL
         CMATDA(ID,IS) = 1./P3
         CMAT6U(ID,IS) =      IMAT6U(ID,IS)*CMATDA(ID,IS)
         CMATUA(ID,IS) = (IMATUA(ID,IS)-P1)*CMATDA(ID,IS)
         BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1

            P1 = ALFA*CMATUA(ID ,ISM)
            P2 = ALFA*CMAT6U(IDM,IS )
            CMAT5L(ID,IS) = IMAT5L(ID,IS)/(1.+P1)
            CMATLA(ID,IS) = IMATLA(ID,IS)/(1.+P2)
            P1 = P1*CMAT5L(ID,IS)
            P2 = P2*CMATLA(ID,IS)
            P3 = IMATDA(ID,IS) + P1 + P2
     &          -CMAT5L(ID,IS)*CMAT6U(ID ,ISM)
     &          -CMATLA(ID,IS)*CMATUA(IDM,IS )
     &          +SMALL
            CMATDA(ID,IS) = 1./P3
            CMAT6U(ID,IS) = (IMAT6U(ID,IS)-P2)*CMATDA(ID,IS)
            CMATUA(ID,IS) = (IMATUA(ID,IS)-P1)*CMATDA(ID,IS)
            BNORM = BNORM + IMATRA(ID,IS)*IMATRA(ID,IS)
         END DO
      END DO
      BNORM = SQRT(BNORM)

      EPSLIN = REPS*BNORM
      UEPS   = 1000.*UNDFLW*BNORM
      IF ( EPSLIN.LT.UEPS .AND. BNORM.GT.0. ) THEN
         IF ( IAMOUT.GE.1 ) THEN
            WRITE (PRINTF,'(A)')
     &         ' ++ SWSIP: the required accuracy is too small'
            WRITE (PRINTF,*)
     &         '           required accuracy    = ',EPSLIN
            WRITE (PRINTF,*)
     &         '           appropriate accuracy = ',UEPS
         END IF
         EPSLIN = UEPS
      END IF

!     --- solve the system by forward and backward substitutions
!         in an iterative manner

 10   IF ( ICONV.EQ.0 .AND. IT.LT.MAXIT ) THEN

         IT    = IT + 1
         ICONV = 1
         RNORM = 0.

         IS  = 1
         ISP = IS + 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)
     &                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &                -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = RES(ID,IS)*CMATDA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)
     &                   -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                   -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                   -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                   -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMATLA(ID,IS)*RES(IDM,IS))*
     &                   CMATDA(ID,IS)
         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)
     &                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMATLA(ID,IS)*RES(IDM,IS))*
     &                CMATDA(ID,IS)

         DO IS = 2, ISSTOP-1
            ISM = IS - 1
            ISP = IS + 1

            IDDUM = IDMIN(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
            IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)
     &                   -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                   -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                   -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                   -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &                   -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID,ISM))*
     &                   CMATDA(ID,IS)

            DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               RES(ID,IS) = IMATRA(ID,IS)
     &                      -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                      -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                      -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                      -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                      -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
               RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
               RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)
     &                                  - CMATLA(ID,IS)*RES(IDM,IS ))*
     &                      CMATDA(ID,IS)
            END DO

            IDDUM = IDMAX(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)
     &                   -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                   -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                   -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &                   -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                   -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)
     &                               - CMATLA(ID,IS)*RES(IDM,IS ))*
     &                   CMATDA(ID,IS)

         END DO

         IS  = ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)
     &                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &                -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID,ISM))*
     &                CMATDA(ID,IS)

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = IMATRA(ID,IS)
     &                   -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                   -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                   -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                   -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
            RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)
     &                               - CMATLA(ID,IS)*RES(IDM,IS ))*
     &                   CMATDA(ID,IS)
         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         RES(ID,IS) = IMATRA(ID,IS)
     &                -IMATDA(ID,IS)*AC2(ID ,IS ,KCGRD(1))
     &                -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &                -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &                -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         RNORM = RNORM + RES(ID,IS)*RES(ID,IS)
         RES(ID,IS) = (RES(ID,IS) - CMAT5L(ID,IS)*RES(ID ,ISM)
     &                            - CMATLA(ID,IS)*RES(IDM,IS ))*
     &                CMATDA(ID,IS)

         IF ( RNORM.GT.1.E8 ) THEN
            IT = MAXIT + 1
            ICONV = 0
            GOTO 10
         END IF
         RNORM=SQRT(RNORM)
         IF ( IAMOUT.EQ.3 .AND. IT.EQ.1 ) RNRM0 = RNORM

         IF ( IAMOUT.EQ.2 ) THEN
            WRITE (PRINTF,'(A,I3,A,E12.6)')
     &            ' ++ SWSIP: iter = ',IT,'    res = ',RNORM
         END IF

         IS    = ISSTOP
         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1

         AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)

         DO IDDUM = IDMAX(IS)-1, IDMIN(IS), -1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            RES(ID,IS) = RES(ID,IS) - CMATUA(ID,IS)*RES(IDP,IS)
            AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)
         END DO

         DO IS = ISSTOP-1, 1, -1
            ISP = IS + 1

            IDDUM = IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1

            RES(ID,IS) = RES(ID,IS) - CMAT6U(ID,IS)*RES(ID,ISP)
            AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)

            DO IDDUM = IDMAX(IS)-1, IDMIN(IS), -1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               RES(ID,IS) = RES(ID,IS) - CMAT6U(ID,IS)*RES(ID ,ISP)
     &                                 - CMATUA(ID,IS)*RES(IDP,IS )
               AC2(ID,IS,KCGRD(1)) = AC2(ID,IS,KCGRD(1)) + RES(ID,IS)
            END DO
         END DO

         IF ( RNORM.GT.UNDFLW**2 .AND. RNORM.GT.EPSLIN ) ICONV = 0
         GOTO 10

      END IF

!     --- investigate the reason to stop

      IF ( ICONV.EQ.0 ) THEN
         AC2(:,:,KCGRD(1)) = AC2OLD(:,:)
         INOCNV = INOCNV + 1
      END IF
      IF ( ICONV.EQ.0 .AND. IAMOUT.GE.1 ) THEN
         IF (ERRPTS.GT.0.AND.IAMMASTER) THEN
            WRITE(ERRPTS,100) IXCGRD(1)+MXF-1, IYCGRD(1)+MYF-1, 2
         END IF
 100     FORMAT (I4,1X,I4,1X,I2)
         WRITE (PRINTF,'(A,I5,A,I5,A)')
     &     ' ++ SWSIP: no convergence in grid point (',
     &      IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
         WRITE (PRINTF,'(A,I3)')
     &     '           total number of iterations     = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '           2-norm of the residual         = ',RNORM
         WRITE (PRINTF,'(A,E12.6)')
     &     '           required accuracy              = ',EPSLIN
      ELSE IF ( IAMOUT.EQ.3 ) THEN
         WRITE (PRINTF,'(A,E12.6)')
     &     ' ++ SWSIP: 2-norm of the initial residual = ',RNRM0
         WRITE (PRINTF,'(A,I3)')
     &     '           total number of iterations     = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '           2-norm of the residual         = ',RNORM
      END IF

!     --- test output

      IF ( TESTFL .AND. ITEST.GE.120 ) THEN
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) '  Subroutine SWSIP'
         WRITE(PRTEST,*)
         WRITE(PRTEST,200) KCGRD(1), MDC, MSC
 200     FORMAT(' SWSIP : POINT MDC MSC              :',3I5)
         WRITE(PRTEST,250) IDDLOW, IDDTOP, ISSTOP
 250     FORMAT(' SWSIP : IDDLOW IDDTOP ISSTOP       :',3I4)
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) ' coefficients of matrix and rhs  '
         WRITE(PRTEST,*)
         WRITE(PRTEST,'(A111)')
     &          ' IS ID         IMATLA         IMATDA'//
     &          '         IMATUA         IMATRA         IMAT5L'//
     &          '         IMAT6U            AC2'
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            DO IS = 1, ISSTOP
               WRITE(PRTEST,300) IS, ID,
     &                           IMATLA(ID,IS), IMATDA(ID,IS),
     &                           IMATUA(ID,IS), IMATRA(ID,IS),
     &                           IMAT5L(ID,IS), IMAT6U(ID,IS),
     &                           AC2(ID,IS,KCGRD(1))
 300           FORMAT(2I3,7E15.7)
            END DO
         END DO
         WRITE(PRTEST,*)
         WRITE(PRTEST,*)'IS ID      LPER          UPER '
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( ID.EQ.1 .OR. ID.EQ.MDC ) THEN
               DO IS = 1, ISSTOP
                  WRITE(PRTEST,350) IS, ID, LOPERI(IS), UPPERI(IS)
 350              FORMAT(2I3,2E15.7)
               END DO
            END IF
         END DO
      END IF

!     --- set matrix coefficients to zero

      IMATDA = 0.
      IMATRA = 0.

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSOR ( AC2   , IMATDA, IMATRA, IMATLA, IMATUA,
     &                   IMAT5L, IMAT6U, AC2OLD, REPS  , MAXIT ,
     &                   IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP,
     &                   IDCMIN, IDCMAX )
!
!****************************************************************
!
      USE SWCOMM1
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4
      USE M_PARALL

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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Nov. 04: New subroutine
!
!  2. Purpose
!
!     Solves penta-diagonal system of equations in
!     spectral space with point SOR method
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     IAMOUT      control parameter indicating the amount of
!                 output required
!                 0: no output
!                 1: only fatal errors will be printed
!                 2: gives output concerning the iteration process
!                 3: additional information about the iteration
!                    is printed
!     IDCMAX      maximum counter in directional space
!     IDCMIN      minimum counter in directional space
!     IDDLOW      minimum direction that is propagated within a sweep
!     IDDTOP      maximum direction that is propagated within a sweep
!     IMAT5L      coefficients of lower diagonal in sigma-space
!     IMAT6U      coefficients of upper diagonal in sigma-space
!     IMATDA      coefficients of main diagonal
!     IMATLA      coefficients of lower diagonal in theta-space
!     IMATUA      coefficients of upper diagonal in theta-space
!     IMATRA      right-hand side
!     INOCNV      integer indicating number of grid points in which
!                 solver does not converged
!     ISSTOP      maximum frequency counter in a sweep
!     MAXIT       the maximum number of iterations to be performed in
!                 the linear solver
!     REPS        relative accuracy of the final approximation
!
      INTEGER IAMOUT, INOCNV, IDDLOW, IDDTOP, ISSTOP, MAXIT
      INTEGER IDCMIN(MSC), IDCMAX(MSC)
      REAL    REPS
      REAL    AC2(MDC,MSC,MCGRD),
     &        IMATDA(MDC,MSC), IMATRA(MDC,MSC),
     &        IMAT5L(MDC,MSC), IMAT6U(MDC,MSC),
     &        IMATLA(MDC,MSC), IMATUA(MDC,MSC),
     &        AC2OLD(MDC,MSC)
!
!  5. Parameter variables
!
!     OMEG  :     relaxation parameter
!
      REAL    OMEG
      PARAMETER (OMEG=0.8)
!
!  6. Local variables
!
!     AC2I  :     intermediate action density
!     ICONV :     indicator for convergence (1=yes, 0=no)
!     ID    :     loop counter in theta-space
!     IDDL  :     minimum counter in theta-space of modulo MDC
!     IDDT  :     maximum counter in theta-space of modulo MDC
!     IDDUM :     loop counter
!     IDINF :     index of point ID with largest error in solution
!     IDM   :     index of point ID-1
!     IDMAX :     local array of maximum counter in theta-space
!     IDMIN :     local array of minimum counter in theta-space
!     IDP   :     index of point ID+1
!     IENT  :     number of entries
!     INVMDA:     inverse of main diagonal
!     IS    :     loop counter in sigma-space
!     ISINF :     index of point IS with largest error in solution
!     ISM   :     index of point IS-1
!     ISP   :     index of point IS+1
!     IT    :     iteration count
!     LOPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!     RES   :     residual
!     RESM  :     inf-norm of residual vector
!     RESM0 :     inf-norm of initial residual vector
!     UPPERI:     auxiliary vector meant for computation in
!                 periodic theta-space
!
      INTEGER ICONV, ID, IDINF, IDDL, IDDT, IDDUM, IDM, IDP, IENT,
     &        IS, ISINF, ISM, ISP, IT
      INTEGER IDMIN(MSC), IDMAX(MSC)
      REAL    AC2I, RES, RESM, RESM0
      REAL    LOPERI(MSC), UPPERI(MSC), INVMDA(MDC,MSC)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWOMPU  (in SWANCOM1)
!
! 12. Structure
!
!     The system of equations is solved using the SOR technique in
!     pointwise manner
!     Note that with omeg=1, the Gauss-Seidel method is recovered
!
!     Convergence is reached, if the difference between two consecutive
!     iteration levels measured w.r.t. the maximum norm is smaller than
!     given tolerance
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSOR')

!     --- initialize arrays

      LOPERI = 0.
      UPPERI = 0.
      INVMDA = 0.

!     --- in case of periodicity in theta-space, store values
!         of matrix coefficients corresponding to left bottom and
!         right top

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).EQ.1 .AND. IDCMAX(IS).EQ.MDC ) THEN
           UPPERI(IS) = IMATLA(  1,IS)
           LOPERI(IS) = IMATUA(MDC,IS)
         END IF
      END DO

!     --- when no bins fall within the sweep, i.e. SECTOR = 0,
!         reset the bounds of sector as 1..MDC (routine SOLPRE
!         has clear the rows in the matrix that do not belong
!         to the sweep)

      DO IS = 1, ISSTOP
         IF ( IDCMIN(IS).LE.IDCMAX(IS) ) THEN
            IDMIN(IS) = IDCMIN(IS)
            IDMAX(IS) = IDCMAX(IS)
         ELSE
            IDMIN(IS) = 1
            IDMAX(IS) = MDC
         END IF
      END DO

!     --- store inverse of main diagonal

      DO IS = 1, ISSTOP
         DO IDDUM = IDMIN(IS), IDMAX(IS)
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( IMATDA(ID,IS).NE.0. ) THEN
               INVMDA(ID,IS) = 1./IMATDA(ID,IS)
            ELSE
               CALL MSGERR ( 3,
     &                     'Main diagonal of spectral matrix is zero!' )
            END IF
         END DO
      END DO

      IT    = 0
      ICONV = 0

!     --- start iteration process

 10   IF ( ICONV.EQ.0 .AND. IT.LT.MAXIT ) THEN

         IT    = IT + 1
         ICONV = 1
         RESM  = 0.
         IDINF = 0
         ISINF = 0

         IS  = 1
         ISP = IS + 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)
     &         -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &         -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &         -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)
     &            -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &            -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &            -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)
     &         -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &         -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &         -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IS = 2, ISSTOP-1
            ISM = IS - 1
            ISP = IS + 1

            IDDUM = IDMIN(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
            IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)
     &            -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &            -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &            -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &            -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

            DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
               ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
               IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
               IDP = MOD ( IDDUM     + MDC , MDC ) + 1

               AC2I = IMATRA(ID,IS)
     &               -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &               -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &               -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &               -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
               AC2I = AC2I*OMEG*INVMDA(ID,IS)+
     &                                     (1.-OMEG)*AC2(ID,IS,KCGRD(1))

               RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
               IF ( RES.GT.RESM ) THEN
                  RESM  = RES
                  IDINF = ID
                  ISINF = IS
               END IF
               AC2(ID,IS,KCGRD(1)) = AC2I

            END DO

            IDDUM = IDMAX(IS)
            ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)
     &            -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &            -IMAT6U(ID,IS)*AC2(ID ,ISP,KCGRD(1))
     &            -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &            -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IS  = ISSTOP
         ISM = IS - 1

         IDDUM = IDMIN(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDP   = MOD ( IDDUM     + MDC , MDC ) + 1
         IDDT  = MOD ( IDMAX(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)
     &         -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &         -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
     &         -UPPERI(IS)*AC2(IDDT,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         DO IDDUM = IDMIN(IS)+1, IDMAX(IS)-1
            ID  = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IDM = MOD ( IDDUM - 2 + MDC , MDC ) + 1
            IDP = MOD ( IDDUM     + MDC , MDC ) + 1

            AC2I = IMATRA(ID,IS)
     &            -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &            -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &            -IMATUA(ID,IS)*AC2(IDP,IS ,KCGRD(1))
            AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

            RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IDINF = ID
               ISINF = IS
            END IF
            AC2(ID,IS,KCGRD(1)) = AC2I

         END DO

         IDDUM = IDMAX(IS)
         ID    = MOD ( IDDUM - 1 + MDC , MDC ) + 1
         IDM   = MOD ( IDDUM - 2 + MDC , MDC ) + 1
         IDDL  = MOD ( IDMIN(IS) - 1 + MDC , MDC ) + 1

         AC2I = IMATRA(ID,IS)
     &         -IMAT5L(ID,IS)*AC2(ID ,ISM,KCGRD(1))
     &         -IMATLA(ID,IS)*AC2(IDM,IS ,KCGRD(1))
     &         -LOPERI(IS)*AC2(IDDL,IS,KCGRD(1))
         AC2I = AC2I*OMEG*INVMDA(ID,IS)+(1.-OMEG)*AC2(ID,IS,KCGRD(1))

         RES = ABS(AC2(ID,IS,KCGRD(1)) - AC2I)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IDINF = ID
            ISINF = IS
         END IF
         AC2(ID,IS,KCGRD(1)) = AC2I

         IF ( RESM.GT.1.E8 ) THEN
            IT = MAXIT + 1
            ICONV = 0
            GOTO 10
         END IF
         IF ( IAMOUT.EQ.2 ) THEN
            WRITE (PRINTF,'(A,I3,A,E12.6,A,I3,A,I3,A)')
     &            ' ++ SWSOR: iter = ',IT,'    res = ',RESM,
     &            ' in (ID,IS) = (',IDINF,',',ISINF,')'
         END IF
         IF ( IAMOUT.EQ.3 .AND. IT.EQ.1 ) RESM0 = RESM

         IF ( RESM.GT.REPS ) ICONV = 0
         GOTO 10

      END IF

!     --- investigate the reason to stop

      IF ( ICONV.EQ.0 ) THEN
         AC2(:,:,KCGRD(1)) = AC2OLD(:,:)
         INOCNV = INOCNV + 1
      END IF
      IF ( ICONV.EQ.0 .AND. IAMOUT.GE.1 ) THEN
         IF (ERRPTS.GT.0.AND.IAMMASTER) THEN
            WRITE(ERRPTS,100) IXCGRD(1)+MXF-1, IYCGRD(1)+MYF-1, 2
         END IF
 100     FORMAT (I4,1X,I4,1X,I2)
         WRITE (PRINTF,'(A,I5,A,I5,A)')
     &     ' ++ SWSOR: no convergence in grid point (',
     &      IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
         WRITE (PRINTF,'(A,I3)')
     &     '           total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '           inf-norm of the residual         = ',RESM
         WRITE (PRINTF,'(A,E12.6)')
     &     '           required accuracy                = ',REPS
      ELSE IF ( IAMOUT.EQ.3 ) THEN
         WRITE (PRINTF,'(A,E12.6)')
     &     ' ++ SWSOR: inf-norm of the initial residual = ',RESM0
         WRITE (PRINTF,'(A,I3)')
     &     '           total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '           inf-norm of the residual         = ',RESM
      END IF

!     --- test output

      IF ( TESTFL .AND. ITEST.GE.120 ) THEN
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) '  Subroutine SWSOR'
         WRITE(PRTEST,*)
         WRITE(PRTEST,200) KCGRD(1), MDC, MSC
 200     FORMAT(' SWSOR : POINT MDC MSC              :',3I5)
         WRITE(PRTEST,250) IDDLOW, IDDTOP, ISSTOP
 250     FORMAT(' SWSOR : IDDLOW IDDTOP ISSTOP       :',3I4)
         WRITE(PRTEST,*)
         WRITE(PRTEST,*) ' coefficients of matrix and rhs  '
         WRITE(PRTEST,*)
         WRITE(PRTEST,'(A111)')
     &          ' IS ID         IMATLA         IMATDA'//
     &          '         IMATUA         IMATRA         IMAT5L'//
     &          '         IMAT6U            AC2'
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            DO IS = 1, ISSTOP
               WRITE(PRTEST,300) IS, ID,
     &                           IMATLA(ID,IS), IMATDA(ID,IS),
     &                           IMATUA(ID,IS), IMATRA(ID,IS),
     &                           IMAT5L(ID,IS), IMAT6U(ID,IS),
     &                           AC2(ID,IS,KCGRD(1))
 300           FORMAT(2I3,7E15.7)
            END DO
         END DO
         WRITE(PRTEST,*)
         WRITE(PRTEST,*)'IS ID      LPER          UPER '
         DO IDDUM = IDDLOW, IDDTOP
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( ID.EQ.1 .OR. ID.EQ.MDC ) THEN
               DO IS = 1, ISSTOP
                  WRITE(PRTEST,350) IS, ID, LOPERI(IS), UPPERI(IS)
 350              FORMAT(2I3,2E15.7)
               END DO
            END IF
         END DO
      END IF

!     --- set matrix coefficients to zero

      IMATDA = 0.
      IMATRA = 0.
      IMATLA = 0.
      IMATUA = 0.
      IMAT5L = 0.
      IMAT6U = 0.

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWMTLB ( N1, N2, M1, M2 )
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
!     40.31: Tim Campbell and John Cazes
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Jul. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Given global loop bounds N1 and N2, compute loop bounds
!     M1 and M2 for calling thread
!
!  4. Argument variables
!
!     M1          lower index of thread loop
!     M2          upper index of thread loop
!     N1          lower index of global loop
!     N2          upper index of global loop
!
      INTEGER N1, N2, M1, M2
!
!  6. Local variables
!
!     ID    :     thread number
!     IENT  :     number of entries
!     NCH   :     auxiliary integer
!     NTH   :     number of threads
!
      INTEGER ID, IENT, NTH, NCH
!$     INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!$     EXTERNAL OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!
!  8. Subroutines used
!
!     ---
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
!     Description of the pseudo code
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWMTLB')

      NTH = 1
      ID  = 0
!$     NTH = OMP_GET_NUM_THREADS()
!$     ID  = OMP_GET_THREAD_NUM()
      NCH = (N2-N1+1)/NTH
      M1  = ID*NCH+N1
      M2  = (ID+1)*NCH+N1-1
      IF(ID.EQ.NTH-1) M2 = N2

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSTPC ( HSACC0, HSACC1, HSACC2, SACC0 , SACC1,
     &                    SACC2 , HSDIFC, TMDIFC, DELHS , DELTM,
     &                    DEP2  , ACCUR , I1MYC , I2MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE M_GENARR
      USE M_PARALL

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
!     40.41: Andre van der Westhuysen
!     40.41: Marcel Zijlema
!     40.93: Andre van der Westhuysen
!
!  1. Updates
!
!     40.41, Jun. 04: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.93, Sep. 08: extended with curvature of Tm
!
!  2. Purpose
!
!     Check convergence based on the relative, absolute
!     and curvature values of wave height and period
!
!  4. Argument variables
!
!     ACCUR       indicates percentage of grid points in
!                 which accuracy is reached
!     DELHS       difference in Hs between last 2 iterations
!     DELTM       difference in Tm between last 2 iterations
!     DEP2        depth
!     HSACC0      significant wave height at iter-2
!     HSACC1      significant wave height at iter-1
!     HSACC2      significant wave height at iter
!     HSDIFC      difference of Hs(i) - Hs(i-2) meant for
!                 computation of curvature of Hs
!     I1MYC       lower index for thread loop over y-grid row
!     I2MYC       upper index for thread loop over y-grid row
!     SACC0       mean wave frequency at iter-2
!     SACC1       mean wave frequency at iter-1
!     SACC2       mean wave frequency at iter
!     TMDIFC      difference of Tm(i) - Tm(i-2) meant for
!                 computation of curvature of Tm
!
      INTEGER I1MYC, I2MYC
      REAL    ACCUR
      REAL    DEP2(MCGRD)          ,
     &        HSACC0(MCGRD)        ,
     &        HSACC1(MCGRD)        ,
     &        HSACC2(MCGRD)        ,
     &        SACC0(MCGRD)         ,
     &        SACC1(MCGRD)         ,
     &        SACC2(MCGRD)         ,
     &        DELHS(MCGRD)         ,
     &        DELTM(MCGRD)         ,
     &        HSDIFC(MCGRD)        ,
     &        TMDIFC(MCGRD)
!
!  6. Local variables
!
!     ACS2  :     auxiliary variable
!     ACS3  :     auxiliary variable
!     HSABS :     absolute value of Hs
!     HSCURV:     curvature value of Hs
!     HSDIFO:     previous value of HSDIFC
!     HSREL :     relative value of Hs
!     IACCUR:     indicates number of grid points in which
!                 accuracy is reached
!     IARR  :     auxiliary array meant for global reduction
!     ID    :     counter of direction
!     IENT  :     number of entries
!     II    :     loop variable
!     INDX  :     index for indirect address
!     IS    :     counter of frequency
!     IX    :     loop counter
!     IX1   :     lower index in x-direction
!     IX2   :     upper index in x-direction
!     IY    :     loop counter
!     IY1   :     lower index in y-direction
!     IY2   :     upper index in y-direction
!     LHEAD :     logical indicating to write header
!     TMABS :     absolute value of Tm
!     TMCURV:     curvature value of Tm
!     TMDIFO:     previous value of TMDIFC
!     TMREL :     relative value of Tm
!     TSTFL :     indicates whether grid point is a test point
!     WETGRD:     number of wet grid points
!     XMOM0 :     zeroth moment
!     XMOM1 :     first moment
!
      INTEGER ID, IS, IENT, II, INDX, IX, IY, IX1, IX2, IY1, IY2
      INTEGER IACCUR, WETGRD, IACCURt, WETGRDt, IARR(2)
      REAL    ACS2, ACS3, HSREL ,HSABS, HSCURV, HSDIFO, TMABS,
     &        TMREL, TMCURV, TMDIFO, XMOM0, XMOM1
      LOGICAL LHEAD, TSTFL
!
!  7. Common blocks used
!
      COMMON/SWSTPC_MT_COM/WETGRD,IACCUR
!
!     SWSTPC_MT_COM    place local summed variables WETGRD and IACCUR
!                      in common block so they will be scoped as shared
!
!  8. Subroutines used
!
!     EQREAL           Boolean function which compares two REAL values
!     STRACE           Tracing routine for debugging
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     SWREDUCE         Performs a global reduction
!
      LOGICAL EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP (in SWANCOM1)
!
! 12. Structure
!
!     master thread initialize the shared variables
!     store Hs and Tm as old values and count number of wet grid points
!     compute new values of Hs and Tm
!     calculate a set of accuracy parameters based on relative,
!         absolute and curvature values of Hs, Tm and check accuracy
!     global sum of IACCUR and WETGRD
!     carry out reductions across all nodes
!
! 13. Source text

      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTPC')

!     --- master thread initialize the shared variables
!$OMP MASTER
      WETGRD = 0
      IACCUR = 0
!$OMP END MASTER
!$OMP BARRIER

      IF ( LMXF ) THEN
         IX1 = 1
      ELSE
         IX1 = 1+IHALOX
      END IF
      IF ( LMXL ) THEN
         IX2 = MXC
      ELSE
         IX2 = MXC-IHALOX
      END IF
      IF ( LMYF ) THEN
         IY1 = I1MYC
      ELSE
         IY1 = 1+IHALOY
      END IF
      IF ( LMYL ) THEN
         IY2 = I2MYC
      ELSE
         IY2 = MYC-IHALOY
      END IF

!     --- store Hs and Tm as old values and count number of wet grid points

      WETGRDt = 0
      DO IX = IX1, IX2
         DO IY = IY1, IY2
            INDX = KGRPNT(IX,IY)
            IF ( DEP2(INDX).GT.DEPMIN ) THEN
               HSACC0(INDX) = MAX( 1.E-20 , HSACC1(INDX) )
               HSACC1(INDX) = MAX( 1.E-20 , HSACC2(INDX) )
               SACC0 (INDX) = MAX( 1.E-20 , SACC1 (INDX) )                40.93
               SACC1 (INDX) = MAX( 1.E-20 , SACC2 (INDX) )
               WETGRDt = WETGRDt + 1
            ELSE
               HSACC0(INDX) = 0.
               HSACC1(INDX) = 0.
               SACC0 (INDX) = 0.                                          40.93
               SACC1 (INDX) = 0.
            END IF
         END DO
      END DO

!     --- compute new values of Hs and Tm

      DO IX = IX1, IX2
         DO IY = IY1, IY2
            INDX = KGRPNT(IX,IY)

            IF ( DEP2(INDX).GT.DEPMIN ) THEN

               XMOM0 = 0.
               XMOM1 = 0.
               DO IS = 1, MSC
                  DO ID = 1, MDC
                     ACS2  = SPCSIG(IS)**2 * AC2(ID,IS,INDX)
                     ACS3  = SPCSIG(IS) * ACS2
                     XMOM0 = XMOM0 + ACS2
                     XMOM1 = XMOM1 + ACS3
                  END DO
               END DO
               XMOM0 = XMOM0 * FRINTF * DDIR
               XMOM1 = XMOM1 * FRINTF * DDIR

               IF ( XMOM0.GT.0. ) THEN
                  HSACC2(INDX) = MAX ( 1.E-20 , 4.*SQRT(XMOM0) )
                  SACC2 (INDX) = MAX ( 1.E-20 , (XMOM1/XMOM0) )
               ELSE
                  HSACC2(INDX) = 1.E-20
                  SACC2 (INDX) = 1.E-20
               END IF

            END IF

         END DO
      END DO

      IACCURt = 0

!     --- calculate a set of accuracy parameters based on relative,
!         absolute and curvature values of Hs and check accuracy

      LHEAD=.TRUE.
      DO IX = IX1, IX2
         DO IY = IY1, IY2
            INDX = KGRPNT(IX,IY)

!           --- determine whether the point is a test point

            TSTFL = .FALSE.
            IF (NPTST.GT.0) THEN
               DO 20 II = 1, NPTST
                  IF (IX.NE.XYTST(2*II-1)) GOTO 20
                  IF (IY.NE.XYTST(2*II  )) GOTO 20
                  TSTFL = .TRUE.
  20           CONTINUE
            END IF

            DELHS(INDX) = 0.0
            DELTM(INDX) = 0.0
            IF ( DEP2(INDX).GT.DEPMIN ) THEN

               HSABS = ABS ( HSACC2(INDX) - HSACC1(INDX) )
               HSREL = HSABS / HSACC2(INDX)
               TMABS = ABS ( (PI2/SACC2(INDX)) - (PI2/SACC1(INDX)) )
               TMREL = TMABS / SACC2(INDX)                                40.93

               HSDIFO       = HSDIFC(INDX)
               HSDIFC(INDX) = 0.5*( HSACC2(INDX) - HSACC0(INDX) )
               HSCURV       = ABS(HSDIFC(INDX) - HSDIFO)/HSACC2(INDX)

               TMDIFO       = TMDIFC(INDX)                                40.93
               TMDIFC(INDX) = 0.5*( SACC2(INDX) - SACC0(INDX) )           40.93
               TMCURV       = ABS(TMDIFC(INDX) - TMDIFO)/SACC2(INDX)      40.93

               DELHS(INDX) = HSABS
               IF (EQREAL(SACC1(INDX),1.E-20) .OR.
     &             EQREAL(SACC2(INDX),1.E-20) ) THEN
                  DELTM(INDX) = 0.
               ELSE
                  DELTM(INDX) = TMABS
               END IF

!              --- add gridpoint in which wave height and period have
!                  reached required accuracy

               IF ( ( HSCURV.LE.PNUMS(15) .AND.
     &               (HSREL.LE.PNUMS(1) .OR. HSABS.LE.PNUMS(2)) ) .AND.   40.93
     &              ( TMCURV.LE.PNUMS(16) .AND.                           40.93
     &               (TMREL.LE.PNUMS(1) .OR. TMABS.LE.PNUMS(3)) ) ) THEN  40.93
                  IACCURt = IACCURt + 1
               END IF

               IF (TSTFL) THEN
                  IF (LHEAD) WRITE(PRINTF,501)
                  WRITE(PRINTF,502) IX+MXF-2, IY+MYF-2, HSABS, HSREL,
     &                              HSCURV, TMABS, TMREL, TMCURV          40.93
 501              FORMAT(25X,'dHabs          ','dHrel          ',
     &                       'Curvature H    ',                           40.93
     &                       'dTabs          ','dTrel          ',         40.93
     &                       'Curvature T    ')                           40.93
 502              FORMAT(1X,SS,'(IX,IY)=(',I5,',',I5,')','  ',
     &                   1PE13.6E2,'  ',1PE13.6E2,'  ',1PE13.6E2,'  ',    40.93
     &                   1PE13.6E2,'  ',1PE13.6E2,'  ',1PE13.6E2)
                  LHEAD=.FALSE.
               END IF

            END IF

         END DO
      END DO
!
!     --- global sum of IACCUR and WETGRD
!$OMP ATOMIC
      IACCUR = IACCUR + IACCURt
!$OMP ATOMIC
      WETGRD = WETGRD + WETGRDt

!     --- carry out reductions across all nodes

!$OMP BARRIER
!$OMP MASTER
      IARR(1) = IACCUR
      IARR(2) = WETGRD
      CALL SWREDUCE ( IARR, 2, SWINT, SWSUM )
      IF (STPNOW()) RETURN
      IACCUR = IARR(1)
      WETGRD = IARR(2)
      ACCUR  = REAL(IACCUR) * 100. / REAL(WETGRD)
!$OMP END MASTER
!$OMP BARRIER
!
!     --- test output
!
!$OMP MASTER
      IF ( ITEST.GE.30 ) THEN
        WRITE(PRINTF,1002) PNUMS(1), PNUMS(2), PNUMS(15)
 1002   FORMAT(' SWSTPC: DHREL DHABS CURV      :',3E12.4)
        WRITE(PRINTF,1008) WETGRD,IACCUR,ACCUR
 1008   FORMAT(' SWSTPC: WETGRD IACCUR ACCUR   :',2I8,E12.4)
      END IF
!$OMP END MASTER

      RETURN
      END
!********************************************************************
!                                                                   *
      SUBROUTINE SETUPP (KGRPNT, MSTPDA, SETPDA, AC2, DEP2, DEPSAV,
     &                   SETUP2, XCGRID, YCGRID, SPCSIG, SPCDIR )         40.41 30.82
!                                                                   *
!********************************************************************
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
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.81: Annette Kieftenburg
!     30.82: IJsbrand Haagsma
!     31.03: Annette Kieftenburg
!     31.04: Nico Booij
!     32.01: Roeland Ris
!     32.03: IJsbrand Haagsma
!     34.01: Jeroen Adema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     32.01, Sept 97: New Subroutine
!     32.03, Feb. 98: Comma added in FORMAT to prevent compilation error
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: transformation of radiation stress in 1D case
!     30.82, Oct. 98: Updated description of several variables
!     30.81, Dec. 98: Argument list KSCIP1 adjusted
!     34.01, Feb. 99: Introducing STPNOW
!     30.82, July 99: Corrected argumentlist SETUPP and SETUP2D
!     30.82, July 99: Corrected argumentlist KSCIP1
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: this routine is reconsidered, cleaned up and moved to SWANCOM1
!
!  2. Purpose
!
!     computes the wave-induced forces and adds the set-up to the depth
!
!  3. Method
!
!     The wave-induced setup is calculated for the one-dimensional
!     mode of SWAN using the following equation:
!
!        d Sxx                d eta
!        ----- +  ( d + eta ) ----- = 0
!         d x                  d x
!
!     This equation is integrated using the forward Euler technique
!
!     For the two-dimensional case, a 2D Poisson equation in general coordinates
!     is solved by means of vertex-centered finite volume method
!
!  4. Argument variables
!
!     AC2       input    action density
!     DEPSAV    input    depth following from original bottom and water level
!     DEP2      i/o      total depth including set-up
!                        on entry: includes previous estimate of set-up
!                        on exit : includes new estimate of set-up
!     KGRPNT    input    indirect addresses for grid points
!     MSTPDA    input    number of (aux.) data per grid point
!                        value is set at 23 in SWANCOM1                   40.41
!     SETPDA    i/o      auxiliary data for computation of set-up
!                        1: x-comp of force, 2: y-comp of force,
!                        3: radiation stress component RSxx, 4: RSxy,
!                        5: RSyy
!                        SETPDA(*,6..MSTPDA) is used as work array
!     SETUP2    output   computed set-up
!     SPCDIR    input    (*,1); spectral directions (radians)             30.82
!                        (*,2); cosine of spectral directions             30.82
!                        (*,3); sine of spectral directions               30.82
!                        (*,4); cosine^2 of spectral directions           30.82
!                        (*,5); cosine*sine of spectral directions        30.82
!                        (*,6); sine^2 of spectral directions             30.82
!     SPCSIG    input    Relative frequencies in computational domain     30.72
!                        in sigma-space                                   30.72
!     XCGRID    input    Coordinates of computational grid in x-direction 30.72
!     YCGRID    input    Coordinates of computational grid in y-direction 30.72
!
      INTEGER MSTPDA, KGRPNT(MXC,MYC)
!
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    DEPSAV(MCGRD)
      REAL    SETPDA(MCGRD,MSTPDA)                                        40.41
      REAL    SETUP2(MCGRD)
      REAL    SPCDIR(MDC,6)                                               30.82
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)                            30.72
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CG          group velocity
!     CHARS       array to pass character info to MSGERR
!     CK          CGO*KWAVE
!     DDET        determinant
!     DEPMAX      maximum depth
!     DEPLOC      local depth
!     DIX         di/dx
!     DIY         di/dy
!     DJX         dj/dx
!     DJY         dj/dy
!     DP1         depth in point i in 1-D case
!     DP2         depth in point i+1 in 1-D case
!     DS2         square of mesh length in x- or y-direction
!     DXI         dx/di
!     DXJ         dx/dj
!     DYI         dy/di
!     DYJ         dy/dj
!     ELOC        local energy
!     ETA1        setup in  point i in 1-D case
!     ETA2        setup in  point i+1 in 1-D case
!     ID          counter in directional space
!     IDXMAX      index of location with maximum depth
!     IENT        number of entries
!     IF1         first non-character in string1
!     IL1         last non-character in string1
!     INDX        address of current grid point (ix,iy)
!     INDXB       address of grid point (ix,iy-1)
!     INDXL       address of grid point (ix-1,iy)
!     INDXR       address of grid point (ix+1,iy)
!     INDXU       address of grid point (ix,iy+1)
!     IS          counter in frequency space
!     IX          counter in x-direction
!     IXLO        counter in x-direction for neighbouring grid point
!     IXUP        counter in x-direction for neighbouring grid point
!     IY          counter in y-direction
!     IYLO        counter in y-direction for neighbouring grid point
!     IYUP        counter in y-direction for neighbouring grid point
!     K           wavenumber
!     LINK        counter for neighbouring grid points
!     MSGSTR      string to pass message to call MSGERR
!     N           CGroup/CPhase
!     ND          derivative of N with respect to depth
!     NEIGHB      boolean variable indicating whether neighbouring point is wet
!     RRDI        1/number of steps in i-direction
!     RRDJ        1/number of steps in j-direction
!     RSXX        xx-component of the radiation stress
!     RSXXI       derivative of RSXX in i-direction
!     RSXXJ       derivative of RSXX in j-direction
!     RSXY        xy-component of the radiation stress
!     RSXYI       derivative of RSXY in i-direction
!     RSXYJ       derivative of RSXY in j-direction
!     RSYY        yy-component of the radiation stress
!     RSYYI       derivative of RSYY in i-direction
!     RSYYJ       derivative of RSYY in j-direction
!     S_UPCOR     total correction to setup (user defined and S_UPDP)
!     S_UPDP      setup at location with maximum depth, before correction
!     SIG         dummy variable for frequency
!     SXX1        radiation stress in  point i in 1-D case
!     SXX2        radiation stress in  point i+1 in 1-D case
!
      INTEGER  IDXMAX, INDXR, INDXU, INDXB,
     &         ID, IENT, INDX, INDXL, IS, IX,
     &         IXLO, IXUP, IY, IYLO, IYUP, LINK
!
      REAL     CK, DDET, DEPLOC, DEPMAX, DIX, DIY, DJX, DJY,              31.04
     &         DP1, DP2, DS2, DXI, DXJ, DYI, DYJ, ELOC, ETA1, ETA2,
     &         RRDI, RRDJ, RSXX,  RSXXI, RSXXJ, RSXY,
     &         RSXYI, RSXYJ, RSYY, RSYYI, RSYYJ,
     &         S_UPCOR,                                                   30.82
     &         S_UPDP,                                                    31.03
     &         SXX1  ,SXX2
      REAL     CG(1), K(1), N(1), ND(1), SIG(1)                           30.82
      INTEGER      IF1, IL1
      CHARACTER*20 INTSTR, CHARS(1)
      CHARACTER*80 MSGSTR
!
      LOGICAL  NEIGHB
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     KSCIP1           Calculates KWAVE, CGO
!     MSGERR           Writes error message
!     SETUP2D          Computation of the change of waterlevel by waves,
!                      a 2D Poisson equation in general coordinates is solved
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOMP
!
! 10. Error messages
!
!     setup in dry point is unequal to zero
!
! 11. Remarks
!
! 12. Structure
!
!     ---------------------------------------------------------
!     For all grid points do
!         If depth > DEPMIN
!         Then Integrate over spectrum to compute RSxx, RSxy, RSyy
!     ---------------------------------------------------------
!     If one-dimensional mode of SWAN
!     Then Calculate Setup in all grid points
!     Else Call SETUP2 to compute setup in all grid points
!     ---------------------------------------------------------
!     S_UPDP is setup in deepest point
!     Add user defined correction to setup                                30.82
!     For all grid points do
!         If dep2 > DEPMIN
!            SETUP2 := SETUP2 - S_UPCOR
!         If dep2 > DEPMIN
!            compute new value for DEP2
!     ---------------------------------------------------------
!     For all grid points do
!         If depth < DEPMIN
!         Then If water level + setup in neighbouring point above
!                   bottom level in current point
!              Then make depth equal to neighbouring water level
!                   + SETUP - bottom level in current point
!     ---------------------------------------------------------
!
! 13. Source text
!
!***********************************************************************

      SAVE     IENT
      DATA     IENT /0/
      CALL STRACE (IENT, 'SETUPP')
!
      DEPMAX = 0.                                                         31.03
      IDXMAX = 0                                                          40.41
!                                                                         31.03
!     --- initializing SETPDA array                                       31.03
!                                                                         31.03
      SETPDA = 0.                                                         40.41
!
      DO IY = 1, MYC
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF (INDX.GT.1) THEN
            IF (DEP2(INDX).GT.DEPMIN) THEN
!                                                                         31.03
!             --- seek deepest point                                      31.03
!                                                                         31.03
              IF (DEPSAV(INDX).GT.DEPMAX) THEN                            31.03
                DEPMAX = DEPSAV(INDX)                                     31.03
                IDXMAX = INDX                                             40.41
              ENDIF                                                       31.03
!
!             --- compute radiation stress components RSXX, RSXY and RSYY
!
              RSXX = 0.
              RSXY = 0.
              RSYY = 0.
              DEPLOC = DEP2(INDX)                                         40.41
              DO IS = 1, MSC
                SIG(1) = SPCSIG(IS)                                       30.82
                CALL KSCIP1 (1,SIG,DEPLOC,K,CG,N,ND)                      30.82
                CK = CG(1) * K(1)                                         30.82
                DO ID = 1, MDC
                  ELOC = SIG(1) * AC2(ID,IS,INDX)                         30.82
!                                  -                                      31.03
!                                  |{cos(Theta)}^2         for i = 4      31.03
!                 SPCDIR(ID,i) is <| sin(Theta)cos(Theta)  for i = 5      31.03
!                                  |{sin(Theta)}^2         for i = 6      31.03
!                                  -                                      31.03
!                                                                         31.03
                  RSXX = RSXX + (CK*SPCDIR(ID,4)+CK - SIG(1)/2.) * ELOC   30.82
                  RSXY = RSXY + CK*SPCDIR(ID,5) * ELOC                    31.03
                  RSYY = RSYY + (CK*SPCDIR(ID,6)+CK - SIG(1)/2.) * ELOC   30.82
                ENDDO
              ENDDO
!
!             --- store radiation stress components in array SETPDA
!
!             DDIR   is width of directional band
!             FRINTF is frequency integration factor df/f
!
              IF (ONED) THEN                                              30.70
!               transform to computational direction
                SETPDA(INDX,3) = DDIR * FRINTF *                          31.04
     &                           ((COSPC*RSXX + SINPC*RSXY) * COSPC +     30.70
     &                            (COSPC*RSXY + SINPC*RSYY) * SINPC)      30.70
              ELSE                                                        30.70
                SETPDA(INDX,3) = RSXX * DDIR * FRINTF
                SETPDA(INDX,4) = RSXY * DDIR * FRINTF
                SETPDA(INDX,5) = RSYY * DDIR * FRINTF
              ENDIF                                                       30.70
            ENDIF                                                         31.03
          ENDIF
        ENDDO
      ENDDO
!
      IF ( ONED ) THEN
!
!       *** compute on the basis of the radiation stresses the setup ***
!
        DO IY = 1, MYC
!          *** boundary condition ***
           SETUP2(KGRPNT(1,IY)) = 0.
           ETA2 = 0.
           DO IX = 1, MXC-1
              INDX  = KGRPNT(IX  ,IY)
              INDXR = KGRPNT(IX+1,IY)
              DP1   = DEP2(INDX )
              DP2   = DEP2(INDXR)
              IF ( INDX .GT.1 .AND. DP1.GT.DEPMIN .AND.                   40.41
     &             INDXR.GT.1 .AND. DP2.GT.DEPMIN ) THEN                  40.41
                 ETA1 = SETUP2(INDX)
                 SXX1 = SETPDA(INDX ,3)                                   40.41 31.04
                 SXX2 = SETPDA(INDXR,3)                                   40.41 31.04
                 ETA2 = ETA1 + ( SXX1 - SXX2 ) / ( 0.5 * ( DP2 + DP1 ) )
                 SETUP2(INDXR) = ETA2
              ELSE
                 SETUP2(INDXR) = 0.                                       40.41
              END IF
           END DO
        END DO
!
      ELSE
!
!       --- compute forces by taking derivative of radiation stress
!
        DO IY = 1, MYC
           DO IX = 1, MXC
              INDX = KGRPNT(IX,IY)                                        31.03
              DEPLOC = DEP2(INDX)
              IF (INDX.GT.1 .AND. DEPLOC.GT.DEPMIN) THEN
                 IF (IX.EQ.1) THEN
                    IXLO = 1
                    IXUP = 2
                 ELSE IF (IX.EQ.MXC) THEN
                    IXLO = MXC-1
                    IXUP = MXC
                 ELSE
                    IXLO = IX-1
                    IXUP = IX+1
                 END IF
                 IF (DEP2(KGRPNT(IXLO,IY)).LE.DEPMIN) IXLO = IX           31.03
                 IF (DEP2(KGRPNT(IXUP,IY)).LE.DEPMIN) IXUP = IX           31.03
                 INDXL = KGRPNT(IXLO,IY)
                 INDXR = KGRPNT(IXUP,IY)
                 IF (IXLO.EQ.IXUP) THEN
                    RRDI = 1.E-20
                 ELSE
                    RRDI = 1. / REAL(IXUP-IXLO)
                 END IF
                 IF (IY.EQ.1) THEN                                        31.03
                    IYLO = 1
                    IYUP = 2
                 ELSE IF (IY.EQ.MYC) THEN
                    IYLO = MYC-1
                    IYUP = MYC
                 ELSE
                    IYLO = IY-1
                    IYUP = IY+1
                 ENDIF
                 IF (DEP2(KGRPNT(IX,IYLO)).LE.DEPMIN) IYLO = IY           31.03
                 IF (DEP2(KGRPNT(IX,IYUP)).LE.DEPMIN) IYUP = IY           31.03
                 INDXB = KGRPNT(IX,IYLO)
                 INDXU = KGRPNT(IX,IYUP)
                 IF (IYLO.EQ.IYUP) THEN
                    RRDJ = 1.E-20
                 ELSE
                    RRDJ = 1. / REAL(IYUP-IYLO)
                 END IF
!
!                --- determine (x,y) derivatives w.r.t. i and j
!
                 DXI = RRDI * (XCGRID(IXUP,IY)-XCGRID(IXLO,IY))
                 DYI = RRDI * (YCGRID(IXUP,IY)-YCGRID(IXLO,IY))
                 DXJ = RRDJ * (XCGRID(IX,IYUP)-XCGRID(IX,IYLO))
                 DYJ = RRDJ * (YCGRID(IX,IYUP)-YCGRID(IX,IYLO))
!
                 RSXXI = RRDI * (SETPDA(INDXR,3)-SETPDA(INDXL,3))
                 RSXXJ = RRDJ * (SETPDA(INDXU,3)-SETPDA(INDXB,3))         31.03
                 RSXYI = RRDI * (SETPDA(INDXR,4)-SETPDA(INDXL,4))
                 RSXYJ = RRDJ * (SETPDA(INDXU,4)-SETPDA(INDXB,4))         31.03
                 RSYYI = RRDI * (SETPDA(INDXR,5)-SETPDA(INDXL,5))
                 RSYYJ = RRDJ * (SETPDA(INDXU,5)-SETPDA(INDXB,5))         31.03
!
                 IF (IXLO.EQ.IXUP.AND.IYLO.EQ.IYUP) THEN                  31.03
!                point surrounded by dry points                           31.03
                    DIX = 0.                                              31.04
                    DIY = 0.                                              31.04
                    DJX = 0.                                              31.04
                    DJY = 0.                                              31.04
                 ELSE IF (IXLO.EQ.IXUP) THEN                              31.03
!                no forces in i-direction                                 31.03
                    DS2 = DXJ**2 + DYJ**2                                 31.04
                    DIX = 0.                                              31.04
                    DIY = 0.                                              31.04
                    DJX = DXJ/DS2                                         31.04
                    DJY = DYJ/DS2                                         31.04
                 ELSE IF (IYLO.EQ.IYUP) THEN                              31.03
!                no forces in j-direction                                 31.03
                    DS2 = DXI**2 + DYI**2                                 31.04
                    DIX = DXI/DS2                                         31.04
                    DIY = DYI/DS2                                         31.04
                    DJX = 0.                                              31.04
                    DJY = 0.                                              31.04
                 ELSE                                                     31.03
!                coefficients for transformation from
!                (i,j)-gradients to (x,y)-gradients
                    DDET = DXI*DYJ - DXJ*DYI
                    DIX  =  DYJ / DDET
                    DIY  = -DXJ / DDET
                    DJX  = -DYI / DDET
                    DJY  =  DXI / DDET
                 END IF                                                   31.04
!
!                --- forces based on spatial gradients of radiation stresses
                 SETPDA(INDX,1) =
     &                  -(RSXXI*DIX + RSXXJ*DJX + RSXYI*DIY + RSXYJ*DJY)  31.03
                 SETPDA(INDX,2) =
     &                  -(RSXYI*DIX + RSXYJ*DJX + RSYYI*DIY + RSYYJ*DJY)  31.03
!
              END IF
           END DO
        END DO
!
!       --- compute set-up in two dimensions
!
         CALL SETUP2D( SETUP2, XCGRID, YCGRID, SETPDA(1,1), SETPDA(1,2),  40.41
     &                 KGRPNT, DEP2, SETPDA(1,6), SETPDA(1,15),           40.41
     &                 SETPDA(1,16) )                                     40.41
!
      END IF
!
      IF (LSETUP.EQ.1) THEN                                               31.04
!       set set-up to 0 for deepest point (This is allowed because the    31.03
!       solution of a Poisson equation + constant is again a solution of  31.03
!       the same Poisson equation)                                        31.03
        S_UPDP = SETUP2(IDXMAX)                                           31.03
        S_UPCOR = S_UPDP - PSETUP(2)                                      30.82
        DO IY = 1, MYC                                                    31.03
           DO IX = 1, MXC                                                 31.03
              INDX = KGRPNT(IX,IY)                                        31.03
              IF (INDX.GT.1) THEN                                         31.03
                 IF (DEP2(INDX).GT.DEPMIN) THEN                           31.03
                    SETUP2(INDX) = SETUP2(INDX) - S_UPCOR                 30.82
                 ELSE
                    IF (ABS(SETUP2(INDX)).GT.1.E-7) THEN                  31.03
                       CHARS(1) = INTSTR(INDX)                            40.41
                       CALL TXPBLA(CHARS(1),IF1,IL1)                      40.41
                       MSGSTR = 'Set-up in dry point with index '//       40.41
     &                          CHARS(1)(IF1:IL1)                         40.41
                       CALL MSGERR ( 2, MSGSTR )                          40.41
                    END IF                                                31.03
                 END IF                                                   31.03
              END IF                                                      31.03
           END DO                                                         31.03
        END DO                                                            31.03
      END IF                                                              31.04
!
!     --- include computed set-up to depth
!
      DO IY = 1, MYC
         DO IX = 1, MXC
            INDX = KGRPNT(IX,IY)
            IF (INDX.GT.1) THEN
               DEP2(INDX) = DEPSAV(INDX) + SETUP2(INDX)
            END IF
         END DO
      END DO
!
!     --- check whether dry points should be inundated
!
      DO IY = 1, MYC
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
!         Note:    KGRPNT(.,.) = 1 means a permanently dry point!         31.03
          IF (INDX.GT.1) THEN
            IF (DEP2(INDX).LE.DEPMIN) THEN
              DO LINK = 1, 4
                NEIGHB = .TRUE.
                IF (LINK.EQ.1) THEN
                  IF (IX.EQ.1) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX-1,IY)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.2) THEN
                  IF (IY.EQ.1) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX,IY-1)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.3) THEN
                  IF (IX.EQ.MXC) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX+1,IY)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ELSE IF (LINK.EQ.4) THEN
                  IF (IY.EQ.MYC) THEN
                    NEIGHB = .FALSE.
                  ELSE
                    INDXL = KGRPNT(IX,IY+1)
                    IF (INDXL.LE.1) NEIGHB = .FALSE.
                  ENDIF
                ENDIF
                IF (NEIGHB) THEN                                          31.03
                  IF (DEPSAV(INDX) + SETUP2(INDXL) .GT. DEPMIN) THEN      31.03
                    SETUP2(INDX) = SETUP2(INDXL)                          31.03
                    DEP2(INDX) = DEPSAV(INDX) + SETUP2(INDXL)             31.03
                  ENDIF                                                   31.03
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
      RETURN
!     end of subroutine SETUPP
      END
!****************************************************************
!
      SUBROUTINE SETUP2D ( SETUP , XCGRID, YCGRID, WFRCX, WFRCY,
     &                     KGRPNT, DEPTH , AMAT  , RHS  , JCTA )
!
!****************************************************************
!
      USE OCPCOMM4
      USE SWCOMM3
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
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Dec. 04: New subroutine
!
!  2. Purpose
!
!     Computation of the change of waterlevel due to waves
!
!  3. Method
!
!     A 2D Poisson equation in general coordinates is solved
!     Vertex-centered finite volume method is employed
!     The system of equations is solved using a SOR method
!
!  4. Argument variables
!
!     AMAT        the coefficient matrix used in the linear system
!     DEPTH       water depth
!     JCTA        Jacobian times contravariant base vectors:
!                                          (K)
!                 JCTA(I,J,K,L)  contains a    in cell with index point I
!                                          L
!                 and point type J
!     KGRPNT      indirect addressing for grid points
!     RHS         the right-hand side vector of the system of equations
!     SETUP       set-up
!     WFRCX       x-component of wave-induced force
!     WFRCY       y-component of wave-induced force
!     XCGRID      x-coordinates of computational grid
!     YCGRID      y-coordinates of computational grid
!
      INTEGER KGRPNT(MXC,MYC)
      REAL    XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      REAL    DEPTH(MCGRD), SETUP(MCGRD), WFRCX(MCGRD), WFRCY(MCGRD)
      REAL    AMAT(MCGRD,9), JCTA(MCGRD,2,2,2), RHS(MCGRD)
!
!  5. Parameter variables
!
!     RELAX :     relaxation parameter
!                 =  -1: relaxation parameter based on gridsizes
!                 <> -1: a fixed (initial) relaxation parameter
!
      REAL    RELAX
      PARAMETER (RELAX=-1.)
!
!  6. Local variables
!
!     CONTRB:     auxiliary variable containing contribution to the matrix
!     DEPF  :     water depth in flux point
!     FACT  :     a factor
!     IAMOUT:     control parameter indicating the amount of
!                 output required
!                 0: no output
!                 1: only fatal errors will be printed
!                 2: gives output concerning the iteration process
!                 3: additional information about the iteration
!                    is printed
!     ICONV :     indicator for convergence (1=yes, 0=no)
!     IENT  :     number of entries
!     II    :     iteration count in case of omega#1
!     INDX  :     index counter for point (ix  ,iy  ) in computational grid
!     INDXB :     index counter for point (ix  ,iy-1) in computational grid
!     INDXL :     index counter for point (ix-1,iy  ) in computational grid
!     INDXLB:     index counter for point (ix-1,iy-1) in computational grid
!     INDXLU:     index counter for point (ix-1,iy+1) in computational grid
!     INDXR :     index counter for point (ix+1,iy  ) in computational grid
!     INDXRB:     index counter for point (ix+1,iy-1) in computational grid
!     INDXRU:     index counter for point (ix+1,iy+1) in computational grid
!     INDXU :     index counter for point (ix  ,iy+1) in computational grid
!     IT    :     iteration count
!     IX    :     counter in x-direction
!     IXINF :     point in x-direction with largest error in solution
!     IY    :     counter in y-direction
!     IYINF :     point in y-direction with largest error in solution
!     JAC   :     Jacobian
!     MAXIT :     the maximum number of iterations to be performed in
!                 the linear solver
!     RES   :     residual
!     RESM  :     inf-norm of residual vector
!     RESMI :     intermediate inf-norm of residual vector
!     RESMO :     inf-norm of residual vector at previous iteration
!     RESM0 :     inf-norm of initial residual vector
!     REPS  :     relative accuracy of the final approximation
!     RHOV  :     estimated largest real eigenvalue
!     SETPI :     intermediate solution for set-up
!     XL    :     measure of convergence speed
!     XOM   :     actual relaxation parameter
!     XOMEG :     optimal overrelaxation parameter for SOR method
!
      INTEGER IAMOUT, ICONV, IENT, INDX, INDXB, INDXL, INDXLB, INDXLU,
     &        INDXR, INDXRB, INDXRU, INDXU, IT, IX, IXINF, IY, IYINF,
     &        II, MAXIT
      REAL    CONTRB, DEPF, FACT, JAC, RES, RESM, RESMI, RESMO, RESM0,
     &        REPS, RHOV, SETPI, XL, XOM, XOMEG
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
!     SETUPP
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     1) point type (=1,2) defines the position of the point in a cell:
!
!        *-------*
!        |       |
!        2       |
!        |       |
!        *---1---*
!
!        * = cell corner, depth point
!
!     2) Neumann boundary condition is imposed on all the boundaries
!
!     3) The determination of the overrelaxation factor by means of
!        alternatively switching between 1 and optimal omega is based
!        on the method as described in
!
!        E.F.F. Botta and M.H.M. Ellenbroek
!        A modified SOR method for the Poisson equation in unsteady
!        free-surface flow calculations
!        J. Comput. Phys., vol. 60, 119-134, 1985
!
! 12. Structure
!
!     initialize some arrays
!     determine contravariant base vectors times Jacobian
!     build right-hand side of the system of equations
!     build the matrix of the linear system
!     in case of nesting put Dirichlet boundary condition
!     set parameters for the solver
!     determine relaxation factor
!     solve the system of equations
!     investigate the reason to stop
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SETUP2D')

!     --- initialize some arrays

      AMAT = 0.
      JCTA = 0.
      RHS  = 0.

!     --- determine contravariant base vector a^(1) times Jacobian
!         in point type 1

      DO IX = 1, MXC-1
         INDX = KGRPNT(IX,1)
         IF ( INDX.GT.1 ) THEN
            JCTA(INDX,1,1,1) = 0.5*( YCGRID(IX+1,2) +
     &                               YCGRID(IX  ,2) -
     &                               YCGRID(IX+1,1) -
     &                               YCGRID(IX  ,1) )
            JCTA(INDX,1,1,2) = 0.5*( XCGRID(IX+1,1) +
     &                               XCGRID(IX  ,1) -
     &                               XCGRID(IX+1,2) -
     &                               XCGRID(IX  ,2) )
         END IF
      END DO
      DO IY = 2, MYC-1
         DO IX = 1, MXC-1
            INDX = KGRPNT(IX,IY)
            IF ( INDX.GT.1 ) THEN
               JCTA(INDX,1,1,1) = 0.25*( YCGRID(IX+1,IY+1) +
     &                                   YCGRID(IX  ,IY+1) -
     &                                   YCGRID(IX+1,IY-1) -
     &                                   YCGRID(IX  ,IY-1) )
               JCTA(INDX,1,1,2) = 0.25*( XCGRID(IX+1,IY-1) +
     &                                   XCGRID(IX  ,IY-1) -
     &                                   XCGRID(IX+1,IY+1) -
     &                                   XCGRID(IX  ,IY+1) )
            END IF
         END DO
      END DO
      DO IX = 1, MXC-1
         INDX = KGRPNT(IX,MYC)
         IF ( INDX.GT.1 ) THEN
            JCTA(INDX,1,1,1) = 0.5*( YCGRID(IX+1,MYC  ) +
     &                               YCGRID(IX  ,MYC  ) -
     &                               YCGRID(IX+1,MYC-1) -
     &                               YCGRID(IX  ,MYC-1) )
            JCTA(INDX,1,1,2) = 0.5*( XCGRID(IX+1,MYC-1) +
     &                               XCGRID(IX  ,MYC-1) -
     &                               XCGRID(IX+1,MYC  ) -
     &                               XCGRID(IX  ,MYC  ) )
         END IF
      END DO

!     --- determine contravariant base vector a^(1) times Jacobian
!         in point type 2

      DO IY = 1, MYC-1
         DO IX = 1, MXC
            INDX = KGRPNT(IX,IY)
            IF ( INDX.GT.1 ) THEN
               JCTA(INDX,2,1,1) = YCGRID(IX,IY+1) - YCGRID(IX,IY  )
               JCTA(INDX,2,1,2) = XCGRID(IX,IY  ) - XCGRID(IX,IY+1)
            END IF
         END DO
      END DO

!     --- determine contravariant base vector a^(2) times Jacobian
!         in point type 1

      DO IY = 1, MYC
         DO IX = 1, MXC-1
            INDX = KGRPNT(IX,IY)
            IF ( INDX.GT.1 ) THEN
               JCTA(INDX,1,2,1) = YCGRID(IX  ,IY) - YCGRID(IX+1,IY)
               JCTA(INDX,1,2,2) = XCGRID(IX+1,IY) - XCGRID(IX  ,IY)
            END IF
         END DO
      END DO

!     --- determine contravariant base vector a^(2) times Jacobian
!         in point type 2

      DO IY = 1, MYC-1
         INDX = KGRPNT(1,IY)
         IF ( INDX.GT.1 ) THEN
            JCTA(INDX,2,2,1) = 0.5*( YCGRID(1,IY+1) +
     &                               YCGRID(1,IY  ) -
     &                               YCGRID(2,IY+1) -
     &                               YCGRID(2,IY  ) )
            JCTA(INDX,2,2,2) = 0.5*( XCGRID(2,IY+1) +
     &                               XCGRID(2,IY  ) -
     &                               XCGRID(1,IY+1) -
     &                               XCGRID(1,IY  ) )
         END IF
         DO IX = 2, MXC-1
            INDX = KGRPNT(IX,IY)
            IF ( INDX.GT.1 ) THEN
               JCTA(INDX,2,2,1) = 0.25*( YCGRID(IX-1,IY+1) +
     &                                   YCGRID(IX-1,IY  ) -
     &                                   YCGRID(IX+1,IY+1) -
     &                                   YCGRID(IX+1,IY  ) )
               JCTA(INDX,2,2,2) = 0.25*( XCGRID(IX+1,IY+1) +
     &                                   XCGRID(IX+1,IY  ) -
     &                                   XCGRID(IX-1,IY+1) -
     &                                   XCGRID(IX-1,IY  ) )
            END IF
         END DO
         INDX = KGRPNT(MXC,IY)
         IF ( INDX.GT.1 ) THEN
            JCTA(INDX,2,2,1) = 0.5*( YCGRID(MXC-1,IY+1) +
     &                               YCGRID(MXC-1,IY  ) -
     &                               YCGRID(MXC  ,IY+1) -
     &                               YCGRID(MXC  ,IY  ) )
            JCTA(INDX,2,2,2) = 0.5*( XCGRID(MXC  ,IY+1) +
     &                               XCGRID(MXC  ,IY  ) -
     &                               XCGRID(MXC-1,IY+1) -
     &                               XCGRID(MXC-1,IY  ) )
         END IF
      END DO

!     --- build right-hand side of the system of equations

!     --- interior domain
      DO IY = 2, MYC-1
         DO IX = 2, MXC-1
            INDX  = KGRPNT(IX  ,IY  )
            INDXL = KGRPNT(IX-1,IY  )
            INDXR = KGRPNT(IX+1,IY  )
            INDXB = KGRPNT(IX  ,IY-1)
            INDXU = KGRPNT(IX  ,IY+1)
            IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!              --- contribution of right flux
               IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
                  RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,1,1,1)
     &                                       *(WFRCX(INDX)+WFRCX(INDXR))
     &                                  - 0.5*JCTA(INDX,1,1,2)
     &                                       *(WFRCY(INDX)+WFRCY(INDXR))
               END IF
!              --- contribution of left flux
               IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
                  RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXL,1,1,1)
     &                                       *(WFRCX(INDXL)+WFRCX(INDX))
     &                                  + 0.5*JCTA(INDXL,1,1,2)
     &                                       *(WFRCY(INDXL)+WFRCY(INDX))
               END IF
!              --- contribution of upper flux
               IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
                  RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,2,2,1)
     &                                       *(WFRCX(INDX)+WFRCX(INDXU))
     &                                  - 0.5*JCTA(INDX,2,2,2)
     &                                       *(WFRCY(INDX)+WFRCY(INDXU))
               END IF
!              --- contribution of bottom flux
               IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
                  RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXB,2,2,1)
     &                                       *(WFRCX(INDXB)+WFRCX(INDX))
     &                                  + 0.5*JCTA(INDXB,2,2,2)
     &                                       *(WFRCY(INDXB)+WFRCY(INDX))
               END IF
            END IF
         END DO
      END DO

!     --- lower boundary (IY=1)
      DO IX = 2, MXC-1
         INDX  = KGRPNT(IX  ,1)
         INDXL = KGRPNT(IX-1,1)
         INDXR = KGRPNT(IX+1,1)
         INDXU = KGRPNT(IX  ,2)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,1,1,1)
     &                                    *(WFRCX(INDX)+WFRCX(INDXR))
     &                               - 0.5*JCTA(INDX,1,1,2)
     &                                    *(WFRCY(INDX)+WFRCY(INDXR))
            END IF
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXL,1,1,1)
     &                                    *(WFRCX(INDXL)+WFRCX(INDX))
     &                               + 0.5*JCTA(INDXL,1,1,2)
     &                                    *(WFRCY(INDXL)+WFRCY(INDX))
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - JCTA(INDX,2,2,1)
     &                                *(WFRCX(INDX)+WFRCX(INDXU))
     &                               - JCTA(INDX,2,2,2)
     &                                *(WFRCY(INDX)+WFRCY(INDXU))
            END IF
         END IF
      END DO

!     --- upper boundary (IY=MYC)
      DO IX = 2, MXC-1
         INDX  = KGRPNT(IX  ,MYC  )
         INDXL = KGRPNT(IX-1,MYC  )
         INDXR = KGRPNT(IX+1,MYC  )
         INDXB = KGRPNT(IX  ,MYC-1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,1,1,1)
     &                                    *(WFRCX(INDX)+WFRCX(INDXR))
     &                               - 0.5*JCTA(INDX,1,1,2)
     &                                    *(WFRCY(INDX)+WFRCY(INDXR))
            END IF
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXL,1,1,1)
     &                                    *(WFRCX(INDXL)+WFRCX(INDX))
     &                               + 0.5*JCTA(INDXL,1,1,2)
     &                                    *(WFRCY(INDXL)+WFRCY(INDX))
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + JCTA(INDXB,2,2,1)
     &                                *(WFRCX(INDXB)+WFRCX(INDX))
     &                               + JCTA(INDXB,2,2,2)
     &                                *(WFRCY(INDXB)+WFRCY(INDX))
            END IF
         END IF
      END DO

!     --- left boundary (IX=1)
      DO IY = 2, MYC-1
         INDX  = KGRPNT(1,IY  )
         INDXR = KGRPNT(2,IY  )
         INDXB = KGRPNT(1,IY-1)
         INDXU = KGRPNT(1,IY+1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - JCTA(INDX,1,1,1)
     &                                *(WFRCX(INDX)+WFRCX(INDXR))
     &                               - JCTA(INDX,1,1,2)
     &                                *(WFRCY(INDX)+WFRCY(INDXR))
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,2,2,1)
     &                                    *(WFRCX(INDX)+WFRCX(INDXU))
     &                               - 0.5*JCTA(INDX,2,2,2)
     &                                    *(WFRCY(INDX)+WFRCY(INDXU))
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXB,2,2,1)
     &                                    *(WFRCX(INDXB)+WFRCX(INDX))
     &                               + 0.5*JCTA(INDXB,2,2,2)
     &                                    *(WFRCY(INDXB)+WFRCY(INDX))
            END IF
         END IF
      END DO

!     --- right boundary (IX=MXC)
      DO IY = 2, MYC-1
         INDX  = KGRPNT(MXC  ,IY  )
         INDXL = KGRPNT(MXC-1,IY  )
         INDXB = KGRPNT(MXC  ,IY-1)
         INDXU = KGRPNT(MXC  ,IY+1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + JCTA(INDXL,1,1,1)
     &                                *(WFRCX(INDXL)+WFRCX(INDX))
     &                               + JCTA(INDXL,1,1,2)
     &                                *(WFRCY(INDXL)+WFRCY(INDX))
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) - 0.5*JCTA(INDX,2,2,1)
     &                                    *(WFRCX(INDX)+WFRCX(INDXU))
     &                               - 0.5*JCTA(INDX,2,2,2)
     &                                    *(WFRCY(INDX)+WFRCY(INDXU))
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               RHS(INDX) = RHS(INDX) + 0.5*JCTA(INDXB,2,2,1)
     &                                    *(WFRCX(INDXB)+WFRCX(INDX))
     &                               + 0.5*JCTA(INDXB,2,2,2)
     &                                    *(WFRCY(INDXB)+WFRCY(INDX))
            END IF
         END IF
      END DO

!     --- left-lower corner (IX=1,IY=1)
      INDX  = KGRPNT(1,1)
      INDXR = KGRPNT(2,1)
      INDXU = KGRPNT(1,2)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of right flux
         IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) - JCTA(INDX,1,1,1)
     &                             *(WFRCX(INDX)+WFRCX(INDXR))
     &                            - JCTA(INDX,1,1,2)
     &                             *(WFRCY(INDX)+WFRCY(INDXR))
         END IF
!        --- contribution of upper flux
         IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) - JCTA(INDX,2,2,1)
     &                             *(WFRCX(INDX)+WFRCX(INDXU))
     &                            - JCTA(INDX,2,2,2)
     &                             *(WFRCY(INDX)+WFRCY(INDXU))
         END IF
      END IF

!     --- right-lower corner (IX=MXC,IY=1)
      INDX  = KGRPNT(MXC  ,1)
      INDXL = KGRPNT(MXC-1,1)
      INDXU = KGRPNT(MXC  ,2)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of left flux
         IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) + JCTA(INDXL,1,1,1)
     &                             *(WFRCX(INDXL)+WFRCX(INDX))
     &                            + JCTA(INDXL,1,1,2)
     &                             *(WFRCY(INDXL)+WFRCY(INDX))
         END IF
!        --- contribution of upper flux
         IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) - JCTA(INDX,2,2,1)
     &                             *(WFRCX(INDX)+WFRCX(INDXU))
     &                            - JCTA(INDX,2,2,2)
     &                             *(WFRCY(INDX)+WFRCY(INDXU))
         END IF
      END IF

!     --- left-upper corner (IX=1,IY=MYC)
      INDX  = KGRPNT(1,MYC  )
      INDXR = KGRPNT(2,MYC  )
      INDXB = KGRPNT(1,MYC-1)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of right flux
         IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) - JCTA(INDX,1,1,1)
     &                             *(WFRCX(INDX)+WFRCX(INDXR))
     &                            - JCTA(INDX,1,1,2)
     &                             *(WFRCY(INDX)+WFRCY(INDXR))
         END IF
!        --- contribution of bottom flux
         IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) + JCTA(INDXB,2,2,1)
     &                             *(WFRCX(INDXB)+WFRCX(INDX))
     &                            + JCTA(INDXB,2,2,2)
     &                             *(WFRCY(INDXB)+WFRCY(INDX))
         END IF
      END IF

!     --- right-upper corner (IX=MXC,IY=MYC)
      INDX  = KGRPNT(MXC  ,MYC  )
      INDXL = KGRPNT(MXC-1,MYC  )
      INDXB = KGRPNT(MXC  ,MYC-1)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of left flux
         IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) + JCTA(INDXL,1,1,1)
     &                             *(WFRCX(INDXL)+WFRCX(INDX))
     &                            + JCTA(INDXL,1,1,2)
     &                             *(WFRCY(INDXL)+WFRCY(INDX))
         END IF
!        --- contribution of bottom flux
         IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
            RHS(INDX) = RHS(INDX) + JCTA(INDXB,2,2,1)
     &                             *(WFRCX(INDXB)+WFRCX(INDX))
     &                            + JCTA(INDXB,2,2,2)
     &                             *(WFRCY(INDXB)+WFRCY(INDX))
         END IF
      END IF

!     --- build the matrix of the linear system

!     --- interior domain
      DO IY = 2, MYC-1
         DO IX = 2, MXC-1
            INDX  = KGRPNT(IX  ,IY  )
            INDXL = KGRPNT(IX-1,IY  )
            INDXR = KGRPNT(IX+1,IY  )
            INDXB = KGRPNT(IX  ,IY-1)
            INDXU = KGRPNT(IX  ,IY+1)
            IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!              --- contribution of right flux
               IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
                  DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXR))
                  JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &                   JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
                  IF ( JAC.NE.0. ) THEN
                     FACT = -DEPF/JAC
                  ELSE
                     FACT = 0.
                  END IF
                  CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                           JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
                  AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
                  AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
                  CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                           JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
                  AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
                  AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
                  AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
                  AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
               END IF
!              --- contribution of left flux
               IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
                  DEPF = 0.5*(DEPTH(INDXL)+DEPTH(INDX))
                  JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &                   JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
                  IF ( JAC.NE.0. ) THEN
                     FACT = -DEPF/JAC
                  ELSE
                     FACT = 0.
                  END IF
                  CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                            JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
                  AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
                  AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
                  CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                            JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
                  AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
                  AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
                  AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
                  AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
               END IF
!              --- contribution of upper flux
               IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
                  DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXU))
                  JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &                   JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
                  IF ( JAC.NE.0. ) THEN
                     FACT = -DEPF/JAC
                  ELSE
                     FACT = 0.
                  END IF
                  CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                           JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
                  AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
                  AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
                  CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                           JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
                  AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
                  AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
                  AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
                  AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
               END IF
!              --- contribution of bottom flux
               IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
                  DEPF = 0.5*(DEPTH(INDXB)+DEPTH(INDX))
                  JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &                   JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
                  IF ( JAC.NE.0. ) THEN
                     FACT = -DEPF/JAC
                  ELSE
                     FACT = 0.
                  END IF
                  CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                            JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
                  AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
                  AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
                  CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                            JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
                  AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
                  AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
                  AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
                  AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
               END IF
            END IF
            IF ( AMAT(INDX,1).EQ.0. ) THEN
               AMAT(INDX,:) = 0.0
               AMAT(INDX,1) = 1.0
               RHS (INDX  ) = 0.0
            END IF
         END DO
      END DO

!     --- lower boundary (IY=1)
      DO IX = 2, MXC-1
         INDX  = KGRPNT(IX  ,1)
         INDXL = KGRPNT(IX-1,1)
         INDXR = KGRPNT(IX+1,1)
         INDXU = KGRPNT(IX  ,2)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXR))
               JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &                JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDXL)+DEPTH(INDX))
               JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &                JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               DEPF = DEPTH(INDX)+DEPTH(INDXU)
               JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &                JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
            AMAT(INDX,5) = AMAT(INDX,5) + 2.*AMAT(INDX,2)
            AMAT(INDX,7) = AMAT(INDX,7) - AMAT(INDX,2)
            AMAT(INDX,2) = 0.
            AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,3)
            AMAT(INDX,8) = AMAT(INDX,8) - AMAT(INDX,3)
            AMAT(INDX,3) = 0.
            AMAT(INDX,6) = AMAT(INDX,6) + 2.*AMAT(INDX,4)
            AMAT(INDX,9) = AMAT(INDX,9) - AMAT(INDX,4)
            AMAT(INDX,4) = 0.
         END IF
         IF ( AMAT(INDX,1).EQ.0. ) THEN
            AMAT(INDX,:) = 0.0
            AMAT(INDX,1) = 1.0
            RHS (INDX  ) = 0.0
         END IF
      END DO

!     --- upper boundary (IY=MYC)
      DO IX = 2, MXC-1
         INDX  = KGRPNT(IX  ,MYC  )
         INDXL = KGRPNT(IX-1,MYC  )
         INDXR = KGRPNT(IX+1,MYC  )
         INDXB = KGRPNT(IX  ,MYC-1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXR))
               JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &                JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDXL)+DEPTH(INDX))
               JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &                JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               DEPF = DEPTH(INDXB)+DEPTH(INDX)
               JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &                JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
            END IF
            AMAT(INDX,5) = AMAT(INDX,5) + 2.*AMAT(INDX,7)
            AMAT(INDX,2) = AMAT(INDX,2) - AMAT(INDX,7)
            AMAT(INDX,7) = 0.
            AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,8)
            AMAT(INDX,3) = AMAT(INDX,3) - AMAT(INDX,8)
            AMAT(INDX,8) = 0.
            AMAT(INDX,6) = AMAT(INDX,6) + 2.*AMAT(INDX,9)
            AMAT(INDX,4) = AMAT(INDX,4) - AMAT(INDX,9)
            AMAT(INDX,9) = 0.
         END IF
         IF ( AMAT(INDX,1).EQ.0. ) THEN
            AMAT(INDX,:) = 0.0
            AMAT(INDX,1) = 1.0
            RHS (INDX  ) = 0.0
         END IF
      END DO

!     --- left boundary (IX=1)
      DO IY = 2, MYC-1
         INDX  = KGRPNT(1,IY  )
         INDXR = KGRPNT(2,IY  )
         INDXB = KGRPNT(1,IY-1)
         INDXU = KGRPNT(1,IY+1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of right flux
            IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
               DEPF = DEPTH(INDX)+DEPTH(INDXR)
               JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &                JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
               CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                        JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXU))
               JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &                JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDXB)+DEPTH(INDX))
               JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &                JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
            END IF
            AMAT(INDX,3) = AMAT(INDX,3) + 2.*AMAT(INDX,2)
            AMAT(INDX,4) = AMAT(INDX,4) - AMAT(INDX,2)
            AMAT(INDX,2) = 0.
            AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,5)
            AMAT(INDX,6) = AMAT(INDX,6) - AMAT(INDX,5)
            AMAT(INDX,5) = 0.
            AMAT(INDX,8) = AMAT(INDX,8) + 2.*AMAT(INDX,7)
            AMAT(INDX,9) = AMAT(INDX,9) - AMAT(INDX,7)
            AMAT(INDX,7) = 0.
         END IF
         IF ( AMAT(INDX,1).EQ.0. ) THEN
            AMAT(INDX,:) = 0.0
            AMAT(INDX,1) = 1.0
            RHS (INDX  ) = 0.0
         END IF
      END DO

!     --- right boundary (IX=MXC)
      DO IY = 2, MYC-1
         INDX  = KGRPNT(MXC  ,IY  )
         INDXL = KGRPNT(MXC-1,IY  )
         INDXB = KGRPNT(MXC  ,IY-1)
         INDXU = KGRPNT(MXC  ,IY+1)
         IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!           --- contribution of left flux
            IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
               DEPF = DEPTH(INDXL)+DEPTH(INDX)
               JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &                JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
               CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                         JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
            END IF
!           --- contribution of upper flux
            IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDX)+DEPTH(INDXU))
               JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &                JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
               AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
               CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                        JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
               AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
               AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
            END IF
!           --- contribution of bottom flux
            IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
               DEPF = 0.5*(DEPTH(INDXB)+DEPTH(INDX))
               JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &                JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
               IF ( JAC.NE.0. ) THEN
                  FACT = -DEPF/JAC
               ELSE
                  FACT = 0.
               END IF
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
               AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
               AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
               CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                         JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
               AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
               AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
               AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
               AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
            END IF
            AMAT(INDX,3) = AMAT(INDX,3) + 2.*AMAT(INDX,4)
            AMAT(INDX,2) = AMAT(INDX,2) - AMAT(INDX,4)
            AMAT(INDX,4) = 0.
            AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,6)
            AMAT(INDX,5) = AMAT(INDX,5) - AMAT(INDX,6)
            AMAT(INDX,6) = 0.
            AMAT(INDX,8) = AMAT(INDX,8) + 2.*AMAT(INDX,9)
            AMAT(INDX,7) = AMAT(INDX,7) - AMAT(INDX,9)
            AMAT(INDX,9) = 0.
         END IF
         IF ( AMAT(INDX,1).EQ.0. ) THEN
            AMAT(INDX,:) = 0.0
            AMAT(INDX,1) = 1.0
            RHS (INDX  ) = 0.0
         END IF
      END DO

!     --- left-lower corner (IX=1,IY=1)
      INDX  = KGRPNT(1,1)
      INDXR = KGRPNT(2,1)
      INDXU = KGRPNT(1,2)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of right flux
         IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDX)+DEPTH(INDXR)
            JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &             JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                     JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
            AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
            CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                     JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
            AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
            AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
            AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
         END IF
!        --- contribution of upper flux
         IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDX)+DEPTH(INDXU)
            JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &             JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                     JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
            AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
            CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                     JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
            AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
            AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
            AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
         END IF
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,3)
         AMAT(INDX,8) = AMAT(INDX,8) - AMAT(INDX,3)
         AMAT(INDX,3) = 0.
         AMAT(INDX,6) = AMAT(INDX,6) + 2.*AMAT(INDX,4)
         AMAT(INDX,9) = AMAT(INDX,9) - AMAT(INDX,4)
         AMAT(INDX,4) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,5)
         AMAT(INDX,6) = AMAT(INDX,6) - AMAT(INDX,5)
         AMAT(INDX,5) = 0.
         AMAT(INDX,8) = AMAT(INDX,8) + 2.*AMAT(INDX,7)
         AMAT(INDX,9) = AMAT(INDX,9) - AMAT(INDX,7)
         AMAT(INDX,7) = 0.
      END IF
      IF ( AMAT(INDX,1).EQ.0. ) THEN
         AMAT(INDX,:) = 0.0
         AMAT(INDX,1) = 1.0
         RHS (INDX  ) = 0.0
      END IF

!     --- right-lower corner (IX=MXC,IY=1)
      INDX  = KGRPNT(MXC  ,1)
      INDXL = KGRPNT(MXC-1,1)
      INDXU = KGRPNT(MXC  ,2)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of left flux
         IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDXL)+DEPTH(INDX)
            JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &             JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                      JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
            AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
            AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
            CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                      JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
            AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
            AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
            AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
         END IF
!        --- contribution of upper flux
         IF ( INDXU.GT.1 .AND. DEPTH(INDXU).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDX)+DEPTH(INDXU)
            JAC  = JCTA(INDX,2,1,1)*JCTA(INDX,2,2,2) -
     &             JCTA(INDX,2,1,2)*JCTA(INDX,2,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,2,1)+
     &                     JCTA(INDX,2,2,2)*JCTA(INDX,2,2,2))
            AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + CONTRB
            CONTRB = FACT*(JCTA(INDX,2,2,1)*JCTA(INDX,2,1,1)+
     &                     JCTA(INDX,2,2,2)*JCTA(INDX,2,1,2))
            AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
            AMAT(INDX,7) = AMAT(INDX,7) - 0.25*CONTRB
            AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
         END IF
         AMAT(INDX,5) = AMAT(INDX,5) + 2.*AMAT(INDX,2)
         AMAT(INDX,7) = AMAT(INDX,7) - AMAT(INDX,2)
         AMAT(INDX,2) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,3)
         AMAT(INDX,8) = AMAT(INDX,8) - AMAT(INDX,3)
         AMAT(INDX,3) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,6)
         AMAT(INDX,5) = AMAT(INDX,5) - AMAT(INDX,6)
         AMAT(INDX,6) = 0.
         AMAT(INDX,8) = AMAT(INDX,8) + 2.*AMAT(INDX,9)
         AMAT(INDX,7) = AMAT(INDX,7) - AMAT(INDX,9)
         AMAT(INDX,9) = 0.
      END IF
      IF ( AMAT(INDX,1).EQ.0. ) THEN
         AMAT(INDX,:) = 0.0
         AMAT(INDX,1) = 1.0
         RHS (INDX  ) = 0.0
      END IF

!     --- left-upper corner (IX=1,IY=MYC)
      INDX  = KGRPNT(1,MYC  )
      INDXR = KGRPNT(2,MYC  )
      INDXB = KGRPNT(1,MYC-1)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of right flux
         IF ( INDXR.GT.1 .AND. DEPTH(INDXR).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDX)+DEPTH(INDXR)
            JAC  = JCTA(INDX,1,1,1)*JCTA(INDX,1,2,2) -
     &             JCTA(INDX,1,1,2)*JCTA(INDX,1,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,1,1)+
     &                     JCTA(INDX,1,1,2)*JCTA(INDX,1,1,2))
            AMAT(INDX,1) = AMAT(INDX,1) - CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + CONTRB
            CONTRB = FACT*(JCTA(INDX,1,1,1)*JCTA(INDX,1,2,1)+
     &                     JCTA(INDX,1,1,2)*JCTA(INDX,1,2,2))
            AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
            AMAT(INDX,4) = AMAT(INDX,4) - 0.25*CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
            AMAT(INDX,9) = AMAT(INDX,9) + 0.25*CONTRB
         END IF
!        --- contribution of bottom flux
         IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDXB)+DEPTH(INDX)
            JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &             JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                      JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
            AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
            AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
            CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                      JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
            AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
            AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
            AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
         END IF
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,8)
         AMAT(INDX,3) = AMAT(INDX,3) - AMAT(INDX,8)
         AMAT(INDX,8) = 0.
         AMAT(INDX,6) = AMAT(INDX,6) + 2.*AMAT(INDX,9)
         AMAT(INDX,4) = AMAT(INDX,4) - AMAT(INDX,9)
         AMAT(INDX,9) = 0.
         AMAT(INDX,3) = AMAT(INDX,3) + 2.*AMAT(INDX,2)
         AMAT(INDX,4) = AMAT(INDX,4) - AMAT(INDX,2)
         AMAT(INDX,2) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,5)
         AMAT(INDX,6) = AMAT(INDX,6) - AMAT(INDX,5)
         AMAT(INDX,5) = 0.
      END IF
      IF ( AMAT(INDX,1).EQ.0. ) THEN
         AMAT(INDX,:) = 0.0
         AMAT(INDX,1) = 1.0
         RHS (INDX  ) = 0.0
      END IF

!     --- right-upper corner (IX=MXC,IY=MYC)
      INDX  = KGRPNT(MXC  ,MYC  )
      INDXL = KGRPNT(MXC-1,MYC  )
      INDXB = KGRPNT(MXC  ,MYC-1)
      IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
!        --- contribution of left flux
         IF ( INDXL.GT.1 .AND. DEPTH(INDXL).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDXL)+DEPTH(INDX)
            JAC  = JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,2) -
     &             JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,1,1)+
     &                      JCTA(INDXL,1,1,2)*JCTA(INDXL,1,1,2))
            AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
            AMAT(INDX,5) = AMAT(INDX,5) - CONTRB
            CONTRB = -FACT*(JCTA(INDXL,1,1,1)*JCTA(INDXL,1,2,1)+
     &                      JCTA(INDXL,1,1,2)*JCTA(INDXL,1,2,2))
            AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
            AMAT(INDX,3) = AMAT(INDX,3) - 0.25*CONTRB
            AMAT(INDX,7) = AMAT(INDX,7) + 0.25*CONTRB
            AMAT(INDX,8) = AMAT(INDX,8) + 0.25*CONTRB
         END IF
!        --- contribution of bottom flux
         IF ( INDXB.GT.1 .AND. DEPTH(INDXB).GT.DEPMIN ) THEN
            DEPF = DEPTH(INDXB)+DEPTH(INDX)
            JAC  = JCTA(INDXB,2,1,1)*JCTA(INDXB,2,2,2) -
     &             JCTA(INDXB,2,1,2)*JCTA(INDXB,2,2,1)
            IF ( JAC.NE.0. ) THEN
               FACT = -DEPF/JAC
            ELSE
               FACT = 0.
            END IF
            CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,2,1)+
     &                      JCTA(INDXB,2,2,2)*JCTA(INDXB,2,2,2))
            AMAT(INDX,1) = AMAT(INDX,1) + CONTRB
            AMAT(INDX,3) = AMAT(INDX,3) - CONTRB
            CONTRB = -FACT*(JCTA(INDXB,2,2,1)*JCTA(INDXB,2,1,1)+
     &                      JCTA(INDXB,2,2,2)*JCTA(INDXB,2,1,2))
            AMAT(INDX,2) = AMAT(INDX,2) - 0.25*CONTRB
            AMAT(INDX,4) = AMAT(INDX,4) + 0.25*CONTRB
            AMAT(INDX,5) = AMAT(INDX,5) - 0.25*CONTRB
            AMAT(INDX,6) = AMAT(INDX,6) + 0.25*CONTRB
         END IF
         AMAT(INDX,5) = AMAT(INDX,5) + 2.*AMAT(INDX,7)
         AMAT(INDX,2) = AMAT(INDX,2) - AMAT(INDX,7)
         AMAT(INDX,7) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,8)
         AMAT(INDX,3) = AMAT(INDX,3) - AMAT(INDX,8)
         AMAT(INDX,8) = 0.
         AMAT(INDX,3) = AMAT(INDX,3) + 2.*AMAT(INDX,4)
         AMAT(INDX,2) = AMAT(INDX,2) - AMAT(INDX,4)
         AMAT(INDX,4) = 0.
         AMAT(INDX,1) = AMAT(INDX,1) + 2.*AMAT(INDX,6)
         AMAT(INDX,5) = AMAT(INDX,5) - AMAT(INDX,6)
         AMAT(INDX,6) = 0.
      END IF
      IF ( AMAT(INDX,1).EQ.0. ) THEN
         AMAT(INDX,:) = 0.0
         AMAT(INDX,1) = 1.0
         RHS (INDX  ) = 0.0
      END IF

!     --- in case of nesting put Dirichlet boundary condition

      IF ( LSETUP.EQ.2 ) THEN
!        --- left and right boundaries
         DO IY = 1, MYC
            DO IX = 1, MXC, MXC-1
               INDX  = KGRPNT(IX,IY)
               AMAT(INDX,:) = 0.0
               AMAT(INDX,1) = 1.0
               IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
                  RHS(INDX) = SETUP(INDX)
               ELSE
                  RHS(INDX) = 0.0
               END IF
            END DO
         END DO
!        --- lower and upper boundaries
         DO IY = 1, MYC, MYC-1
            DO IX = 1, MXC
               INDX  = KGRPNT(IX,IY)
               AMAT(INDX,:) = 0.0
               AMAT(INDX,1) = 1.0
               IF ( INDX.GT.1 .AND. DEPTH(INDX).GT.DEPMIN ) THEN
                  RHS(INDX) = SETUP(INDX)
               ELSE
                  RHS(INDX) = 0.0
               END IF
            END DO
         END DO
      END IF

!     --- set parameters for the solver

      REPS   = PNUMS(23)
      IAMOUT = INT(PNUMS(24))
      MAXIT  = INT(PNUMS(25))

      CSETUP = .TRUE.

!     --- determine relaxation factor

      IF ( RELAX.EQ.-1. ) THEN
         IF ( MCGRD.LT.100 ) THEN
            RHOV = 1. - 1./REAL(MCGRD)
         ELSE IF ( MCGRD.LT.1000 ) THEN
            RHOV = 1. - 3./REAL(MCGRD)
         ELSE
            RHOV = 1. - 10./REAL(MCGRD)
         END IF
         XOMEG = 2./(1.+SQRT(1.-RHOV*RHOV))
      ELSE
         XOMEG = RELAX
      END IF

      IT    = 0
      ICONV = 0
      RESM  = 1.
      XOM   = 1.

!     --- solve the system of equations

 10   IF ( ICONV.EQ.0 .AND. IT.LT.MAXIT ) THEN

         IT    = IT + 1
         ICONV = 1

         RESMO = RESM
         RESM  = 0.
         IXINF = 0
         IYINF = 0

!        --- interior domain
         DO IY = 2, MYC-1
            DO IX = 2, MXC-1
               INDX   = KGRPNT(IX  ,IY  )
               INDXL  = KGRPNT(IX-1,IY  )
               INDXR  = KGRPNT(IX+1,IY  )
               INDXB  = KGRPNT(IX  ,IY-1)
               INDXU  = KGRPNT(IX  ,IY+1)
               INDXLB = KGRPNT(IX-1,IY-1)
               INDXRB = KGRPNT(IX+1,IY-1)
               INDXLU = KGRPNT(IX-1,IY+1)
               INDXRU = KGRPNT(IX+1,IY+1)

               SETPI = RHS(INDX) - AMAT(INDX,2)*SETUP(INDXLB)
     &                           - AMAT(INDX,3)*SETUP(INDXB )
     &                           - AMAT(INDX,4)*SETUP(INDXRB)
     &                           - AMAT(INDX,5)*SETUP(INDXL )
     &                           - AMAT(INDX,6)*SETUP(INDXR )
     &                           - AMAT(INDX,7)*SETUP(INDXLU)
     &                           - AMAT(INDX,8)*SETUP(INDXU )
     &                           - AMAT(INDX,9)*SETUP(INDXRU)
               SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

               RES = ABS(SETUP(INDX) - SETPI)
               IF ( RES.GT.RESM ) THEN
                  RESM  = RES
                  IXINF = IX
                  IYINF = IY
               END IF
               SETUP(INDX) = SETPI

            END DO
         END DO

!        --- lower boundary (IY=1)
         DO IX = 2, MXC-1
            INDX   = KGRPNT(IX  ,1)
            INDXL  = KGRPNT(IX-1,1)
            INDXR  = KGRPNT(IX+1,1)
            INDXU  = KGRPNT(IX  ,2)
            INDXLU = KGRPNT(IX-1,2)
            INDXRU = KGRPNT(IX+1,2)

            SETPI = RHS(INDX) - AMAT(INDX,5)*SETUP(INDXL )
     &                        - AMAT(INDX,6)*SETUP(INDXR )
     &                        - AMAT(INDX,7)*SETUP(INDXLU)
     &                        - AMAT(INDX,8)*SETUP(INDXU )
     &                        - AMAT(INDX,9)*SETUP(INDXRU)
            SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

            RES = ABS(SETUP(INDX) - SETPI)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IXINF = IX
               IYINF = 1
            END IF
            SETUP(INDX) = SETPI

         END DO

!        --- upper boundary (IY=MYC)
         DO IX = 2, MXC-1
            INDX   = KGRPNT(IX  ,MYC  )
            INDXL  = KGRPNT(IX-1,MYC  )
            INDXR  = KGRPNT(IX+1,MYC  )
            INDXB  = KGRPNT(IX  ,MYC-1)
            INDXLB = KGRPNT(IX-1,MYC-1)
            INDXRB = KGRPNT(IX+1,MYC-1)

            SETPI = RHS(INDX) - AMAT(INDX,2)*SETUP(INDXLB)
     &                        - AMAT(INDX,3)*SETUP(INDXB )
     &                        - AMAT(INDX,4)*SETUP(INDXRB)
     &                        - AMAT(INDX,5)*SETUP(INDXL )
     &                        - AMAT(INDX,6)*SETUP(INDXR )
            SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

            RES = ABS(SETUP(INDX) - SETPI)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IXINF = IX
               IYINF = MYC
            END IF
            SETUP(INDX) = SETPI

         END DO

!        --- left boundary (IX=1)
         DO IY = 2, MYC-1
            INDX   = KGRPNT(1,IY  )
            INDXR  = KGRPNT(2,IY  )
            INDXB  = KGRPNT(1,IY-1)
            INDXU  = KGRPNT(1,IY+1)
            INDXRB = KGRPNT(2,IY-1)
            INDXRU = KGRPNT(2,IY+1)

            SETPI = RHS(INDX) - AMAT(INDX,3)*SETUP(INDXB )
     &                        - AMAT(INDX,4)*SETUP(INDXRB)
     &                        - AMAT(INDX,6)*SETUP(INDXR )
     &                        - AMAT(INDX,8)*SETUP(INDXU )
     &                        - AMAT(INDX,9)*SETUP(INDXRU)
            SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

            RES = ABS(SETUP(INDX) - SETPI)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IXINF = 1
               IYINF = IY
            END IF
            SETUP(INDX) = SETPI

         END DO

!        --- right boundary (IX=MXC)
         DO IY = 2, MYC-1
            INDX   = KGRPNT(MXC  ,IY  )
            INDXL  = KGRPNT(MXC-1,IY  )
            INDXB  = KGRPNT(MXC  ,IY-1)
            INDXU  = KGRPNT(MXC  ,IY+1)
            INDXLB = KGRPNT(MXC-1,IY-1)
            INDXLU = KGRPNT(MXC-1,IY+1)

            SETPI = RHS(INDX) - AMAT(INDX,2)*SETUP(INDXLB)
     &                        - AMAT(INDX,3)*SETUP(INDXB )
     &                        - AMAT(INDX,5)*SETUP(INDXL )
     &                        - AMAT(INDX,7)*SETUP(INDXLU)
     &                        - AMAT(INDX,8)*SETUP(INDXU )
            SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

            RES = ABS(SETUP(INDX) - SETPI)
            IF ( RES.GT.RESM ) THEN
               RESM  = RES
               IXINF = MXC
               IYINF = IY
            END IF
            SETUP(INDX) = SETPI

         END DO

!        --- left-lower corner (IX=1,IY=1)
         INDX   = KGRPNT(1,1)
         INDXR  = KGRPNT(2,1)
         INDXU  = KGRPNT(1,2)
         INDXRU = KGRPNT(2,2)

         SETPI = RHS(INDX) - AMAT(INDX,6)*SETUP(INDXR )
     &                     - AMAT(INDX,8)*SETUP(INDXU )
     &                     - AMAT(INDX,9)*SETUP(INDXRU)
         SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

         RES = ABS(SETUP(INDX) - SETPI)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IXINF = 1
            IYINF = 1
         END IF
         SETUP(INDX) = SETPI

!        --- right-lower corner (IX=MXC,IY=1)
         INDX   = KGRPNT(MXC  ,1)
         INDXL  = KGRPNT(MXC-1,1)
         INDXU  = KGRPNT(MXC  ,2)
         INDXLU = KGRPNT(MXC-1,2)

         SETPI = RHS(INDX) - AMAT(INDX,5)*SETUP(INDXL )
     &                     - AMAT(INDX,7)*SETUP(INDXLU)
     &                     - AMAT(INDX,8)*SETUP(INDXU )
         SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

         RES = ABS(SETUP(INDX) - SETPI)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IXINF = MXC
            IYINF = 1
         END IF
         SETUP(INDX) = SETPI

!        --- left-upper corner (IX=1,IY=MYC)
         INDX   = KGRPNT(1,MYC  )
         INDXR  = KGRPNT(2,MYC  )
         INDXB  = KGRPNT(1,MYC-1)
         INDXRB = KGRPNT(2,MYC-1)

         SETPI = RHS(INDX) - AMAT(INDX,3)*SETUP(INDXB )
     &                     - AMAT(INDX,4)*SETUP(INDXRB)
     &                     - AMAT(INDX,6)*SETUP(INDXR )
         SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

         RES = ABS(SETUP(INDX) - SETPI)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IXINF = 1
            IYINF = MYC
         END IF
         SETUP(INDX) = SETPI

!        --- right-upper corner (IX=MXC,IY=MYC)
         INDX   = KGRPNT(MXC  ,MYC  )
         INDXL  = KGRPNT(MXC-1,MYC  )
         INDXB  = KGRPNT(MXC  ,MYC-1)
         INDXLB = KGRPNT(MXC-1,MYC-1)

         SETPI = RHS(INDX) - AMAT(INDX,2)*SETUP(INDXLB)
     &                     - AMAT(INDX,3)*SETUP(INDXB )
     &                     - AMAT(INDX,5)*SETUP(INDXL )
         SETPI = XOM*SETPI/AMAT(INDX,1) + (1.-XOM)*SETUP(INDX)

         RES = ABS(SETUP(INDX) - SETPI)
         IF ( RES.GT.RESM ) THEN
            RESM  = RES
            IXINF = MXC
            IYINF = MYC
         END IF
         SETUP(INDX) = SETPI

         IF ( RESM.GT.1.E8 ) THEN
            IT = MAXIT + 1
            ICONV = 0
            GOTO 10
         END IF
         IF ( IAMOUT.EQ.2 ) THEN
            WRITE (PRINTF,'(A,I3,A,E12.6,A,I3,A,I3,A)')
     &            ' ++ SETUP2D: iter = ',IT,'    res = ',RESM,
     &            ' in (IX,IY) = (',IXINF,',',IYINF,')'
         END IF
         IF ( IAMOUT.EQ.3 .AND. IT.EQ.1 ) RESM0 = RESM

         IF ( RESM.GT.REPS ) THEN
            IF ( XOM.NE.1. ) THEN
               II = II + 1
               IF ( II.EQ.0 ) RESMI = RESM
               IF ( II.GT.0 .AND.
     &              RESM.GT.10.*RESMI*(XOMEG-1.)**II ) XOM = 1.
            ELSE
               XL = RESM/RESMO
               IF ( XL.GT.0.9*RHOV*RHOV ) THEN
                  XOM = XOMEG
                  II  = -10
               END IF
            END IF
            ICONV = 0
         END IF

         GOTO 10

      END IF

!     --- investigate the reason to stop

      IF ( ICONV.EQ.0 ) THEN
         CSETUP = .FALSE.
         IF ( RESM.GT.1.E8 ) SETUP = 0.
      END IF
      IF ( ICONV.EQ.0 .AND. IAMOUT.GE.1 ) THEN
         WRITE (PRINTF,'(A)')
     &     ' ++ SETUP2D: no convergence for 2D Poisson equation'
         WRITE (PRINTF,'(A,I4)')
     &     '             total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '             inf-norm of the residual         = ',RESM
         WRITE (PRINTF,'(A,E12.6)')
     &     '             required accuracy                = ',REPS
      ELSE IF ( IAMOUT.EQ.3 ) THEN
         WRITE (PRINTF,'(A,E12.6)')
     &     ' ++ SETUP2D: inf-norm of the initial residual = ',RESM0
         WRITE (PRINTF,'(A,I4)')
     &     '             total number of iterations       = ',IT
         WRITE (PRINTF,'(A,E12.6)')
     &     '             inf-norm of the residual         = ',RESM
      END IF

      RETURN
      END
