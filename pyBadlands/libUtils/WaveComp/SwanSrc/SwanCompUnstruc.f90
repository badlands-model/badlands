subroutine SwanCompUnstruc ( ac2, ac1, compda, spcsig, spcdir, xytst, cross, it )
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
!   Authors
!
!   40.80: Marcel Zijlema
!   40.85: Marcel Zijlema
!   40.95: Marcel Zijlema
!   41.02: Marcel Zijlema
!   41.07: Marcel Zijlema
!   41.10: Marcel Zijlema
!   41.20: Casey Dietrich
!
!   Updates
!
!   40.80,     July 2007: New subroutine
!   40.85,   August 2008: add propagation, generation and redistribution terms for output purposes
!   40.95,     June 2008: parallelization of unSWAN using MESSENGER of ADCIRC
!   41.02, February 2009: implementation of diffraction
!   41.07,   August 2009: bug fix: never-ending sweep is prevented
!   41.10,   August 2009: parallelization using OpenMP directives
!   41.20,     June 2010: extension to tightly coupled ADCIRC+SWAN model
!
!   Purpose
!
!   Performs one time step for solution of wave action equation on unstructured grid
!
!   Method
!
!   A vertex-based algorithm is employed in which the variables are stored at the vertices of the mesh
!   The equation is solved in each vertex assuming a constant spectral grid resolution in all vertices
!   The propagation terms in both geographic and spectral spaces are integrated implicitly
!   Sources are treated explicitly and sinks implicitly
!   The calculation of the source terms is carried out in the original SWAN routines, e.g., SOURCE
!
!   The wave action equation is solved iteratively
!   A number of iterations are carried out until convergence is reached
!   In each iteration, a number of sweeps are carried out
!   Per sweep, a loop over all vertices is executed
!   The solution of each vertex must be updated geographically before proceeding to the next one
!   The two upwave faces connecting the vertex to be updated enclose those wave directions that can be processed in the spectral space
!   A sweep is complete when all vertices are updated geographically (but not necessarily in whole spectral space)
!   The process continues by sweeping through the vertices in a reverse manner, allowing waves to propagate from other directions
!   An iteration is complete when all vertices are updated in both geographic and spectral spaces
!
!   Modules used
!
    use ocpcomm4
    use swcomm1
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
    use m_parall
!PUN    use SIZES, only: SZ, MNPROC
!PUN    use MESSENGER
!ADC    use couple2adcirc, only: MakeBoundariesReflective
!ADC    use NodalAttributes, only: FoundSwanWaveRefrac, LoadSwanWaveRefrac, SwanWaveRefrac
!
    implicit none
!
!   Argument variables
!
    integer, dimension(nfaces), intent(in)         :: cross  ! contains sequence number of obstacles for each face
                                                             ! where they crossing or zero if no crossing
    integer, intent(in)                            :: it     ! counter in time space
    integer, dimension(NPTST), intent(in)          :: xytst  ! test points for output purposes
!
    real, dimension(MDC,MSC,nverts), intent(in)    :: ac1    ! action density at previous time level
    real, dimension(MDC,MSC,nverts), intent(inout) :: ac2    ! action density at current time level
    real, dimension(nverts,MCMVAR), intent(inout)  :: compda ! array containing space-dependent info (e.g. depth)
    real, dimension(MDC,6), intent(in)             :: spcdir ! (*,1): spectral direction bins (radians)
                                                             ! (*,2): cosine of spectral directions
                                                             ! (*,3): sine of spectral directions
                                                             ! (*,4): cosine^2 of spectral directions
                                                             ! (*,5): cosine*sine of spectral directions
                                                             ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)               :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell     ! cell index / loop counter
    integer                               :: id        ! loop counter over direction bins
    integer                               :: iddlow    ! minimum direction bin that is propagated within a sweep
    integer                               :: iddtop    ! maximum direction bin that is propagated within a sweep
    integer, parameter                    :: idebug=0  ! level of debug output:
                                                       ! 0 = no output
                                                       ! 1 = print extra output for debug purposes
    integer                               :: idtot     ! maximum number of bins in directional space for considered sweep
    integer                               :: idwmax    ! maximum counter for spectral wind direction
    integer                               :: idwmin    ! minimum counter for spectral wind direction
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: iface     ! face index
    integer                               :: inocnt    ! inocnv counter for calling thread
    integer                               :: inocnv    ! integer indicating number of vertices in which solver does not converged
    integer                               :: is        ! loop counter over frequency bins
    integer                               :: isslow    ! minimum frequency that is propagated within a sweep
    integer                               :: isstop    ! maximum frequency that is propagated within a sweep
    integer                               :: istat     ! indicate status of allocation
    integer                               :: istot     ! maximum number of bins in frequency space for considered sweep
    integer                               :: iter      ! iteration counter
    integer                               :: ivert     ! vertex index
    integer                               :: ivlow     ! lower index in range of vertices in calling thread
    integer                               :: ivup      ! upper index in range of vertices in calling thread
    integer                               :: j         ! loop counter
    integer                               :: jc        ! loop counter
    integer                               :: k         ! loop counter
    integer                               :: kend      ! last index in range of vertices
    integer                               :: kstart    ! first index in range of vertices
    integer                               :: kstep     ! step through range of vertices depending on sweep direction
    integer                               :: kvert     ! loop counter over vertices
    integer                               :: l         ! counter
    integer, dimension(2)                 :: link      ! local face connected to considered vertex where an obstacle crossed
    integer                               :: mnisl     ! minimum sigma-index occured in applying limiter
    integer                               :: mxnfl     ! maximum number of use of limiter in spectral space
    integer                               :: mxnfr     ! maximum number of use of rescaling in spectral space
    integer                               :: npfl      ! number of vertices in which limiter is used
    integer                               :: npfr      ! number of vertices in which rescaling is used
    integer                               :: nupdv     ! number of updates of present vertex (based on surrounding cells)
    integer                               :: swpcol    ! own sweep color (0/1=black/white) in calling thread
    integer                               :: swpdir    ! sweep counter
    integer                               :: swpnr     ! sweep number
    integer                               :: tid       ! thread number
    integer, dimension(3)                 :: v         ! vertices in present cell
    integer                               :: vb        ! vertex of begin of present face
    integer                               :: ve        ! vertex of end of present face
    integer, dimension(2)                 :: vu        ! upwave vertices in present cell
    integer, dimension(24)                :: wwint     ! counters for 4 wave-wave interactions (see routine FAC4WW)
    !
    integer, dimension(:), allocatable    :: idcmax    ! maximum frequency-dependent counter in directional space
    integer, dimension(:), allocatable    :: idcmin    ! minimum frequency-dependent counter in directional space
    integer, dimension(:), allocatable    :: iscmax    ! maximum direction-dependent counter in frequency space
    integer, dimension(:), allocatable    :: iscmin    ! minimum direction-dependent counter in frequency space
    integer, dimension(:), allocatable    :: islmin    ! lowest sigma-index occured in applying limiter
    integer, dimension(:), allocatable    :: nflim     ! number of frequency use of limiter in each vertex
    integer, dimension(:), allocatable    :: nrscal    ! number of frequency use of rescaling in each vertex
!$  integer, dimension(:), allocatable    :: tlist     ! vertex list for calling thread
    !
    real                                  :: abrbot    ! near bottom excursion
    real                                  :: accur     ! percentage of active vertices in which required accuracy has been reached
    real                                  :: acnrmo    ! norm of difference of previous iteration
    real, dimension(2)                    :: acnrms    ! array containing infinity norms
    real                                  :: dal1      ! a coefficent for the 4 wave-wave interactions
    real                                  :: dal2      ! another coefficent for the 4 wave-wave interactions
    real                                  :: dal3      ! just another coefficent for the 4 wave-wave interactions
!PUN    real(SZ), dimension(1)                :: dum1      ! a dummy real meant for UPDATER
!PUN    real(SZ), dimension(1)                :: dum2      ! a dummy real meant for UPDATER
    real                                  :: dummy     ! dummy variable (to be used in existing SWAN routine call)
    real                                  :: etot      ! total wave energy density
    real                                  :: fpm       ! Pierson Moskowitz frequency
    real                                  :: frac      ! fraction of total active vertices
    real                                  :: hm        ! maximum wave height
    real                                  :: hs        ! significant wave height
    real                                  :: kmespc    ! mean average wavenumber based on the WAM formulation
    real                                  :: nwetp     ! total number of active vertices
    real                                  :: qbloc     ! fraction of breaking waves
    real, dimension(2)                    :: rdx       ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2)                    :: rdy       ! second component of contravariant base vector rdy(b) = a^(b)_2
    real                                  :: rhof      ! asymptotic convergence factor
    real                                  :: rval1     ! a dummy value
    real                                  :: rval2     ! a dummy value
    real                                  :: smebrk    ! mean frequency based on the first order moment
    real                                  :: snlc1     ! a coefficent for the 4 wave-wave interactions
    real                                  :: stopcr    ! stopping criterion for stationary solution
    real                                  :: thetaw    ! mean direction of the wind speed vector with respect to ambient current
    real                                  :: ufric     ! wind friction velocity
    real, dimension(5)                    :: usrset    ! auxiliary array to store user-defined settings of 3rd generation mode
    real                                  :: wind10    ! magnitude of the wind speed vector with respect to ambient current
    real, dimension(8)                    :: wwawg     ! weight coefficients for the 4 wave-wave interactions (see routine FAC4WW)
    real, dimension(8)                    :: wwswg     ! weight coefficients for the 4 wave-wave interactions semi-implicitly (see routine FAC4WW)
    real                                  :: xis       ! difference between succeeding frequencies for computing 4 wave-wave interactions
    !
    real, dimension(:,:), allocatable     :: ac2old    ! array to store action density before solving system of equations
!PUN    real(SZ), dimension(:), allocatable   :: ac2loc    ! local action density meant for UPDATER
    real, dimension(:,:), allocatable     :: alimw     ! maximum energy by wind growth
                                                       ! this auxiliary array is used because the maximum value has to be checked
                                                       ! direct after solving the action balance equation
    real, dimension(:,:,:), allocatable   :: amat      ! coefficient matrix of system of equations in spectral space
    real, dimension(:,:), allocatable     :: rhs       ! right-hand side of system of equations in spectral space
    real, dimension(:,:), allocatable     :: cad       ! wave transport velocity in theta-direction
    real, dimension(:,:), allocatable     :: cas       ! wave transport velocity in sigma-direction
    real, dimension(:,:,:), allocatable   :: cax       ! wave transport velocity in x-direction
    real, dimension(:,:,:), allocatable   :: cay       ! wave transport velocity in y-direction
    real, dimension(:,:), allocatable     :: cgo       ! group velocity
    real, dimension(:,:), allocatable     :: da1c      ! implicit interaction contribution of first quadruplet, current bin (unfolded space)
    real, dimension(:,:), allocatable     :: da1m      ! implicit interaction contribution of first quadruplet, current bin -1 (unfolded space)
    real, dimension(:,:), allocatable     :: da1p      ! implicit interaction contribution of first quadruplet, current bin +1 (unfolded space)
    real, dimension(:,:), allocatable     :: da2c      ! implicit interaction contribution of second quadruplet, current bin (unfolded space)
    real, dimension(:,:), allocatable     :: da2m      ! implicit interaction contribution of second quadruplet, current bin -1 (unfolded space)
    real, dimension(:,:), allocatable     :: da2p      ! implicit interaction contribution of second quadruplet, current bin +1 (unfolded space)
    real, dimension(:,:,:), allocatable   :: disc0     ! explicit part of dissipation in present vertex for output purposes
    real, dimension(:,:,:), allocatable   :: disc1     ! implicit part of dissipation in present vertex for output purposes
    real, dimension(:,:), allocatable     :: dsnl      ! total interaction contribution of quadruplets to the main diagonal matrix
    real, dimension(:,:,:), allocatable   :: genc0     ! explicit part of generation in present vertex for output purposes
    real, dimension(:,:,:), allocatable   :: genc1     ! implicit part of generation in present vertex for output purposes
    real, dimension(:), allocatable       :: hscurr    ! wave height at current iteration level
    real, dimension(:), allocatable       :: hsdifc    ! difference in wave height of current and one before previous iteration
    real, dimension(:), allocatable       :: hsprev    ! wave height at previous iteration level
    real, dimension(:,:), allocatable     :: kwave     ! wave number
    real, dimension(:,:), allocatable     :: leakcf    ! leak coefficient in present vertex for output purposes
    real, dimension(:,:,:), allocatable   :: memnl4    ! auxiliary array to store results of 4 wave-wave interactions in full spectral space
    real, dimension(:,:,:), allocatable   :: obredf    ! action reduction coefficient based on transmission
    real, dimension(:,:,:), allocatable   :: redc0     ! explicit part of redistribution in present vertex for output purposes
    real, dimension(:,:,:), allocatable   :: redc1     ! implicit part of redistribution in present vertex for output purposes
    real, dimension(:,:), allocatable     :: reflso    ! contribution to the source term due to reflection
    real, dimension(:,:), allocatable     :: sa1       ! explicit interaction contribution of first quadruplet (unfolded space)
    real, dimension(:,:), allocatable     :: sa2       ! explicit interaction contribution of second quadruplet (unfolded space)
    real, dimension(:,:), allocatable     :: sfnl      ! total interaction contribution of quadruplets to the right-hand side
    real, dimension(:,:,:,:), allocatable :: swtsda    ! several source terms computed at test points
    real, dimension(:), allocatable       :: tmcurr    ! mean period at current iteration level
    real, dimension(:), allocatable       :: tmdifc    ! difference in mean period of current and one before previous iteration
    real, dimension(:), allocatable       :: tmprev    ! mean period at previous iteration level
    real, dimension(:,:,:), allocatable   :: trac0     ! explicit part of propagation in present vertex for output purposes
    real, dimension(:,:,:), allocatable   :: trac1     ! implicit part of propagation in present vertex for output purposes
    real, dimension(:,:), allocatable     :: ue        ! energy density for computing 4 wave-wave interactions (unfolded space)
    !
    logical                               :: fguess    ! indicate whether first guess need to be applied or not
    logical                               :: lpredt    ! indicate whether action density in first iteration need to be estimated or not
    logical                               :: swpfull   ! indicate whether all necessary sweeps are done or not
    !
    logical, dimension(:,:), allocatable  :: anybin    ! true if bin is active in considered sweep
    logical, dimension(:,:), allocatable  :: anyblk    ! true if bin is blocked by a counter current based on a CFL criterion
    logical, dimension(:), allocatable    :: anywnd    ! true if wind input is active in considered bin
    logical, dimension(:,:), allocatable  :: groww     ! check for each frequency whether the waves are growing or not
                                                       ! in a spectral direction
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes
    !
!$  integer, external                     :: omp_get_num_threads ! number of OpenMP threads being used
!$  integer, external                     :: omp_get_thread_num  ! get thread number
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanCompUnstruc')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! some initializations
    !
    ICMAX  = 3         ! stencil size
    PROPSL = PROPSC
    !
    IXCGRD(1) = -9999  ! to be used in routines SINTGRL and SOURCE so that Ursell number and
    IYCGRD(1) = -9999  ! quadruplets are calculated once in each vertex during an iteration
    !
    tid   = 0
    ivlow = 1
    ivup  = nverts
    !
    ! print all the settings used in SWAN run
    !
    if ( it == 1 .and. ITEST > 0 ) call SWPRSET (spcsig)
    !
    ! print test points
    !
    if ( NPTST > 0 ) then
       do j = 1, NPTST
          write (PRINTF,107) j, xytst(j)
       enddo
    endif
    !
    ! allocation of shared arrays
    !
!TIMG    call SWTSTA(101)
    allocate(islmin(nverts))
    allocate( nflim(nverts))
    allocate(nrscal(nverts))
    !
    allocate(hscurr(nverts))
    allocate(hsprev(nverts))
    allocate(hsdifc(nverts))
    allocate(tmcurr(nverts))
    allocate(tmprev(nverts))
    allocate(tmdifc(nverts))
!PUN    allocate(ac2loc(nverts))
    !
    allocate(swtsda(MDC,MSC,NPTSTA,MTSVAR))
    !
    if ( IQUAD > 2 ) then
       allocate(memnl4(MDC,MSC,nverts), stat = istat)
       if ( istat /= 0 ) then
          call msgerr ( 4, 'Allocation problem in SwanCompUnstruc: array memnl4 ' )
          return
       endif
    else
       allocate(memnl4(0,0,0))
    endif
    !
    !$ allocate(tlist(nverts))
!TIMG    call SWTSTO(101)
    !
    ! initialization of shared arrays
    !
    hscurr = 0.
    hsprev = 0.
    hsdifc = 0.
    tmcurr = 0.
    tmprev = 0.
    tmdifc = 0.
    !
    swtsda = 0.
    !
    ! spawn a parallel region
    !
    !$omp parallel default(shared) &
    !$omp private(cad, cas, cax, cay, cgo, kwave) &
    !$omp private(idcmin, idcmax, iscmin, iscmax, anybin) &
    !$omp private(amat, rhs, ac2old) &
    !$omp private(anywnd, obredf, reflso, alimw, groww, anyblk) &
    !$omp private(disc0, disc1, genc0, genc1, redc0, redc1, trac0, trac1, leakcf) &
    !$omp private(ue, sa1, sa2, sfnl) &
    !$omp private(da1c, da1p, da1m, da2c, da2p, da2m, dsnl) &
    !$omp private(tid, ivlow, ivup, iter, swpcol, kvert, kstart, kend, kstep) &
    !$omp private(ivert, nupdv, jc, k, j, icell, v, vu, l, swpnr, rdx, rdy, lpredt, vb, ve, iface, link, inocnt) &
    !$omp private(iddlow, iddtop, idtot, isslow, isstop, istot) &
    !$omp private(abrbot, kmespc, idwmin, idwmax, hs, etot, qbloc, ufric, fpm, thetaw, hm, wind10, smebrk) &
    !$omp copyin(ICMAX, COSLAT, IPTST, TESTFL)
    !
    ! print number of threads set by environment
    !
    !$omp master
    !$    if ( it == 1 ) &
    !$    write (SCREEN,'(a,i2/)') ' Number of threads during execution of parallel region = ', omp_get_num_threads()
    !$omp end master
    !
    ! get thread number
    !
    !$ tid = omp_get_thread_num()
    tid = tid + 1
    !
    ! allocation of private arrays
    !
!TIMG    call SWTSTA(101)
    allocate(   cad(MDC,MSC      ))
    allocate(   cas(MDC,MSC      ))
    allocate(   cax(MDC,MSC,ICMAX))
    allocate(   cay(MDC,MSC,ICMAX))
    allocate (  cgo(    MSC,ICMAX))
    allocate (kwave(    MSC,ICMAX))
    !
    allocate(idcmax(    MSC))
    allocate(idcmin(    MSC))
    allocate(iscmax(MDC    ))
    allocate(iscmin(MDC    ))
    allocate(anybin(MDC,MSC))
    !
    allocate(  amat(MDC,MSC,5))
    allocate(   rhs(MDC,MSC  ))
    allocate(ac2old(MDC,MSC  ))
    !
    allocate(anywnd(MDC))
    allocate(obredf(MDC,MSC,2))
    allocate(reflso(MDC,MSC))
    allocate( alimw(MDC,MSC))
    allocate( groww(MDC,MSC))
    allocate(anyblk(MDC,MSC))
    !
    allocate( disc0(MDC,MSC,MDISP))
    allocate( disc1(MDC,MSC,MDISP))
    allocate( genc0(MDC,MSC,MGENR))
    allocate( genc1(MDC,MSC,MGENR))
    allocate( redc0(MDC,MSC,MREDS))
    allocate( redc1(MDC,MSC,MREDS))
    allocate( trac0(MDC,MSC,MTRNP))
    allocate( trac1(MDC,MSC,MTRNP))
    allocate(leakcf(MDC,MSC      ))
    !
    ! calculate ranges of spectral space for arrays related to 4 wave-wave interactions
    !
    !$omp single
!TIMG    call SWTSTA(135)
    if ( IQUAD > 0 ) call FAC4WW ( xis, snlc1, dal1, dal2, dal3, spcsig, wwint, wwawg, wwswg )
!TIMG    call SWTSTO(135)
    !$omp end single
    !
    if ( IQUAD > 0 ) then
       allocate(  ue(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate( sa1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate( sa2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate(sfnl(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       if ( IQUAD == 1 ) then
          allocate(da1c(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da1p(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da1m(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2c(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2p(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2m(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(dsnl(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       else
          allocate(da1c(0,0))
          allocate(da1p(0,0))
          allocate(da1m(0,0))
          allocate(da2c(0,0))
          allocate(da2p(0,0))
          allocate(da2m(0,0))
          allocate(dsnl(0,0))
       endif
    else
       allocate(  ue(0,0))
       allocate( sa1(0,0))
       allocate( sa2(0,0))
       allocate(sfnl(0,0))
       allocate(da1c(0,0))
       allocate(da1p(0,0))
       allocate(da1m(0,0))
       allocate(da2c(0,0))
       allocate(da2p(0,0))
       allocate(da2m(0,0))
       allocate(dsnl(0,0))
    endif
!TIMG    call SWTSTO(101)
    !
    ! marks vertices active and non-active
    !
    !$omp single
    nwetp = 0.
    do kvert = 1, nverts
       if ( compda(kvert,JDP2) > DEPMIN ) then
          vert(kvert)%active = .true.
          nwetp = nwetp +1.
       else
          vert(kvert)%active = .false.
       endif
    enddo
    if ( it == 1 .and. ITEST > 0 ) write (PRINTF,108) nint(nwetp), nwetp*100./real(nverts)
    !
    ! First guess of action density will be applied if 3rd generation mode is employed and wind is active (IWIND > 2)
    ! Note: this first guess is not used in nonstationary run (NSTATC > 0) or hotstart (ICOND = 4)
    !
    if ( IWIND > 2 .and. NSTATC == 0 .and. ICOND /= 4 ) then
       fguess = .true.
    else
       fguess = .false.
    endif
    !$omp end single
    !
    ! compute load-balanced spatial loop bounds for each thread
    !
    !$ call SwanThreadBounds( nwetp, ivlow, ivup, tlist )
    !
!TIMG    call SWTSTA(103)
    iterloop: do iter = 1, ITERMX
       !
       ! some initializations
       !
       if ( IQUAD  > 2 ) memnl4(1:MDC,1:MSC,ivlow:ivup) = 0.
       if ( ITRIAD > 0 &
                       ) compda(ivlow:ivup,JURSEL) = 0.
       !
       compda(ivlow:ivup,JDISS) = 0.
       compda(ivlow:ivup,JLEAK) = 0.
       compda(ivlow:ivup,JDSXB) = 0.
       compda(ivlow:ivup,JDSXS) = 0.
       compda(ivlow:ivup,JDSXW) = 0.
       compda(ivlow:ivup,JDSXV) = 0.
       compda(ivlow:ivup,JGENR) = 0.
       compda(ivlow:ivup,JGSXW) = 0.
       compda(ivlow:ivup,JREDS) = 0.
       compda(ivlow:ivup,JRSXQ) = 0.
       compda(ivlow:ivup,JRSXT) = 0.
       compda(ivlow:ivup,JTRAN) = 0.
       compda(ivlow:ivup,JTSXG) = 0.
       compda(ivlow:ivup,JTSXT) = 0.
       compda(ivlow:ivup,JTSXS) = 0.
       compda(ivlow:ivup,JRADS) = 0.
       compda(ivlow:ivup,JQB  ) = 0.
       !
       inocnt = 0
       !
       !$omp master
       inocnv = 0
       !
       islmin = 9999
       nflim  = 0
       nrscal = 0
       !
       acnrms = -9999.
       !
       ! During first iteration, first guess of action density is based on 2nd generation mode
       ! After first iteration, user-defined settings are re-used
       !
       if ( fguess ) then
          !
          if ( iter == 1 )then
             !
             ! save user-defined settings of 3rd generation mode
             ! Note: surf breaking, bottom friction and triads may be still active
             !
             usrset(1) = IWIND
             usrset(2) = IWCAP
             usrset(3) = IQUAD
             usrset(4) = PNUMS(20)
             usrset(5) = PNUMS(30)
             !
             ! first guess settings
             !
             IWIND     = 2               ! if first guess should be based on 1st generation mode, set IWIND = 1
             IWCAP     = 0
             IQUAD     = 0
             PNUMS(20) = 1.E22           ! no limiter
             PNUMS(30) = 0.              ! no under-relaxation
             !
             write (PRINTF,101)
             !
          elseif ( iter == 2 ) then
             !
             ! re-activate user-defined settings of 3rd generation mode
             !
             IWIND     = usrset(1)
             IWCAP     = usrset(2)
             IQUAD     = usrset(3)
             PNUMS(20) = usrset(4)
             PNUMS(30) = usrset(5)
             !
             write (PRINTF,102)
             !
          endif
          !
          ! print info
          !
          if ( iter < 3 ) then
             write (PRINTF,103) iter, PNUMS(20), PNUMS(30)
             write (PRINTF,104) IWIND, IWCAP, IQUAD
             write (PRINTF,105) ISURF, IBOT , ITRIAD
             write (PRINTF,106) IVEG
          endif
          !
       endif
       !
       ! calculate diffraction parameter and its derivatives
       !
       if ( IDIFFR /= 0 ) call SwanDiffPar ( ac2, compda(1,JDP2), spcsig )
       !
       ! all vertices are set untagged except non-active ones
       !
       vert(:)%fullupdated = .false.
       !
       do kvert = 1, nverts
          do jc = 1, vert(kvert)%noc         ! all cells around vertex
             vert(kvert)%updated(jc) = 0
          enddo
       enddo
       !
       do kvert = 1, nverts
          !
          ! in case of non-active vertex set action density equal to zero
          !
          if ( .not.vert(kvert)%active ) then
             !
             vert(kvert)%fullupdated = .true.
             !
             do jc = 1, vert(kvert)%noc
                icell = vert(kvert)%cell(jc)%atti(CELLID)
                vert(kvert)%updated(jc) = icell
             enddo
             !
             ac2(:,:,kvert) = 0.
             !
          endif
          !
       enddo
       !
       swpdir = 0
       !$omp end master
       !
       ! synchronize threads before loop over sweep directions
       !$omp barrier
       !
       ! loop over a number of sweeps through grid
       !
       sweeploop: do
          !
          ! check whether all vertices in grid are completely updated (in both geographic and spectral spaces)
          !
          !$omp single
          swpfull = .true.
          do kvert = 1, nverts
             if ( .not.vert(kvert)%fullupdated ) then
                swpfull = .false.
                exit
             endif
          enddo
          !$omp end single
          !
          if ( swpfull ) exit sweeploop
          !
          !$omp single
          swpdir = swpdir + 1
          !
          if ( SCREEN /= PRINTF ) then
             if ( NSTATC == 1 ) then
                if ( IAMMASTER ) write (SCREEN,110) CHTIME, it, iter, swpdir
             else
                write (PRINTF,120) iter, swpdir
                if ( IAMMASTER ) then
                   if ( swpdir == 1 ) then
                      write (SCREEN,120) iter, swpdir
                   else
                      write (SCREEN,130) iter, swpdir
                   endif
                endif
             endif
          endif
          !$omp end single
          !
          ! loop over vertices in the grid (depending on sweep direction)
          !
          swpcol = mod(swpdir,2)
          !
          if ( mod(tid,2) == swpcol ) then
             kstart = ivlow
             kend   = ivup
             kstep  = 1
          else
             kstart = ivup
             kend   = ivlow
             kstep  = -1
          endif
          !
          vertloop: do kvert = kstart, kend, kstep
            !
            ivert = vlist(kvert)
            !$ ivert = tlist(kvert)
!ADC            !
!ADC            ! allow SWAN to handle wave refraction as a nodal attribute
!ADC            if ( LoadSwanWaveRefrac .and. FoundSwanWaveRefrac ) then
!ADC               IREFR = nint(SwanWaveRefrac(ivert))
!ADC            endif
            !
            if ( vert(ivert)%active ) then   ! this active vertex needs to be updated
               !
               ! determine whether the present vertex is a test point
               !
               IPTST  = 0
               TESTFL = .false.
               if ( NPTST > 0 ) then
                  do j = 1, NPTST
                     if ( ivert /= xytst(j) ) cycle
                     IPTST  = j
                     TESTFL = .true.
                  enddo
               endif
               !
               nupdv = 0
               !
               celloop: do jc = 1, vert(ivert)%noc
                  !
                  icell = vert(ivert)%cell(jc)%atti(CELLID)
                  !
                  if ( vert(ivert)%updated(jc) == icell ) then   ! this vertex in present cell is already updated
                     nupdv = nupdv + 1
                     cycle celloop
                  endif
                  !
                  v(1) = cell(icell)%atti(CELLV1)
                  v(2) = cell(icell)%atti(CELLV2)
                  v(3) = cell(icell)%atti(CELLV3)
                  !
                  ! pick up two upwave vertices
                  !
                  do k = 1, 3
                     if ( v(k) == ivert ) then
                        vu(1) = v(mod(k  ,3)+1)
                        vu(2) = v(mod(k+1,3)+1)
                        exit
                     endif
                  enddo
                  !
                  l = 0
                  do k = 1, 2
                     do j = 1, vert(vu(k))%noc
                        if ( vert(vu(k))%updated(j) /= 0 .or. vert(vu(k))%atti(VMARKER) /= 0 ) then   ! this upwave vertex is geographically updated or known
                           l = l + 1
                           exit
                        endif
                     enddo
                  enddo
                  !
                  if ( l /= 2 ) cycle celloop   ! the present vertex cannot be updated
                  !
                  ! stores vertices of computational stencil
                  !
                  vs(1) = ivert
                  vs(2) = vu(1)
                  vs(3) = vu(2)
                  !
                  KCGRD = vs    ! to be used in some original SWAN routines
                  !
                  swpnr = 0                                              ! this trick assures to calculate Ursell number and
                  if ( all(mask=vert(ivert)%updated(:)==0) ) swpnr = 1   ! quadruplets only once in each vertex during an iteration
                  !
                  ! compute wavenumber and group velocity in points of stencil
                  !
!TIMG                  call SWTSTA(110)
                  call SwanDispParm ( kwave, cgo, compda(1,JDP2), spcsig )
!TIMG                  call SWTSTO(110)
                  !
                  ! compute wave transport velocities in points of stencil for all directions
                  !
!TIMG                  call SWTSTA(111)
                  call SwanPropvelX ( cax, cay, compda(1,JVX2), compda(1,JVY2), cgo, spcdir(1,2), spcdir(1,3) )
!TIMG                  call SWTSTO(111)
                  !
                  ! compute local contravariant base vectors at present vertex
                  !
                  do k = 1, 3
                     if ( v(k) == ivert ) then
                        rdx(1) = cell(icell)%geom(k)%rdx1
                        rdx(2) = cell(icell)%geom(k)%rdx2
                        rdy(1) = cell(icell)%geom(k)%rdy1
                        rdy(2) = cell(icell)%geom(k)%rdy2
                        exit
                     endif
                  enddo
                  !
                  ! in case of spherical coordinates determine cosine of latitude (in degrees)
                  ! and recalculate local contravariant base vectors
                  !
                  if ( KSPHER > 0 ) then
                     do k = 1, ICMAX
                        COSLAT(k) = cos(DEGRAD*(vert(vs(k))%attr(VERTY) + YOFFS))
                     enddo
                     do j = 1, 2
                        rdx(j) = rdx(j) / (COSLAT(1) * LENDEG)
                        rdy(j) = rdy(j) / LENDEG
                     enddo
                  endif
                  !
                  ! compute spectral directions for the considered sweep in present vertex
                  !
!TIMG                  call SWTSTA(112)
                  call SwanSweepSel ( idcmin, idcmax, anybin, iscmin, iscmax, &
                                      iddlow, iddtop, idtot , isslow, isstop, &
                                      istot , cax   , cay   , rdx   , rdy   , &
                                      spcsig)
!TIMG                  call SWTSTO(112)
                  !
                  if ( idtot > 0 ) then
                     !
                     ! compute propagation velocities in spectral space for the considered sweep in present vertex
                     !
!TIMG                     call SWTSTA(113)
                     call SwanPropvelS ( cad           , cas           , compda(1,JVX2), compda(1,JVY2), &
                                         compda(1,JDP1), compda(1,JDP2), cax           , cay           , &
                                         kwave         , cgo           , spcsig        , idcmin        , &
                                         idcmax        , spcdir(1,2)   , spcdir(1,3)   , spcdir(1,4)   , &
                                         spcdir(1,5)   , spcdir(1,6)   , rdx           , rdy           , &
                                         jc            )
!TIMG                     call SWTSTO(113)
                     !
                     ! estimate action density in case of first iteration at cold start in stationary mode
                     ! (since it is zero in first stationary run)
                     !
                     lpredt = .false.
                     if ( iter == 1 .and. ICOND /= 4 .and. NSTATC == 0 ) then
                        lpredt = .true.
                        compda(ivert,JHS) = 0.
                        goto 20
                     endif
 10                  if ( lpredt ) then
                        call SPREDT (swpdir, ac2   , cax, cay, idcmin, idcmax, &
                                     isstop, anybin, rdx, rdy, obredf)
                        lpredt = .false.
                     endif
                     !
                     ! calculate various integral parameters
                     !
!TIMG                     call SWTSTA(116)
                     call SINTGRL (spcdir, kwave        , ac2  , compda(1,JDP2), qbloc , compda(1,JURSEL), &
                                   rdx   , rdy          , dummy, etot          , abrbot, compda(1,JUBOT) , &
                                   hs    , compda(1,JQB), hm   , kmespc        , smebrk, compda(1,JPBOT) , &
                                   swpnr )
!TIMG                     call SWTSTO(116)
                     !
                     compda(ivert,JHS) = hs
 20                  continue
                     !
                     ! compute transmission and/or reflection if obstacle is present in computational stencil
                     !
                     obredf = 1.
                     reflso = 0.
                     !
!TIMG                     call SWTSTA(136)
                     if ( NUMOBS > 0 ) then
                        !
                        ! determine obstacle for the link(s) in the computational stencil
                        !
                        do j = 1, cell(icell)%nof
                           !
                           ! determine vertices of the local face
                           !
                           vb = cell(icell)%face(j)%atti(FACEV1)
                           ve = cell(icell)%face(j)%atti(FACEV2)
                           !
                           if ( vb==ivert .or. ve==ivert ) then
                              !
                              if ( vb==vu(1) .or. ve==vu(1) ) then
                                 !
                                 iface = cell(icell)%face(j)%atti(FACEID)
                                 link(1) = cross(iface)
                                 !
                              elseif ( vb==vu(2) .or. ve==vu(2) ) then
                                 !
                                 iface = cell(icell)%face(j)%atti(FACEID)
                                 link(2) = cross(iface)
                                 !
                              endif
                              !
                           endif
                           !
                        enddo
                        !
                        if ( link(1)/=0 .or. link(2)/=0 ) then
                           !
                           call SWTRCF ( compda(1,JWLV2), compda(1,JHS), link, obredf, ac2, reflso, dummy , dummy  , &
                                         dummy          , cax          , cay , rdx   , rdy, anybin, spcsig, spcdir )
                           !
                        endif
                        !
                     endif
!TIMG                     call SWTSTO(136)
                     if (lpredt) goto 10
                     !
                     ! compute the transport part of the action balance equation
                     !
!TIMG                     call SWTSTA(118)
                     call SwanTranspAc ( amat  , rhs   , leakcf, ac2   , ac1   , &
                                         cgo   , cax   , cay   , cad   , cas   , &
                                         anybin, rdx   , rdy   , spcsig, spcdir, &
                                         obredf, idcmin, idcmax, iscmin, iscmax, &
                                         iddlow, iddtop, isslow, isstop, anyblk, &
                                         trac0 , trac1 )
!TIMG                     call SWTSTO(118)
                     !
                     ! compute the source part of the action balance equation
                     !
                     if ( .not.OFFSRC ) then
                        !
                        ! initialize wind friction velocity and Pierson Moskowitz frequency
                        !
                        ufric = 1.e-15
                        fpm   = 1.e-15
                        !
                        ! compute the wind speed, mean wind direction, the PM frequency,
                        ! wind friction velocity and the minimum and maximum counters for
                        ! active wind input
                        !
!TIMG                        call SWTSTA(115)
                        if ( IWIND > 0 ) call WINDP1 ( wind10, thetaw, idwmin        , idwmax        , &
                                                       fpm   , ufric , compda(1,JWX2), compda(1,JWY2), &
                                                       anywnd, spcdir, compda(1,JVX2), compda(1,JVY2), &
                                                       spcsig )
!ADC                                                       spcsig, ivert )
!TIMG                        call SWTSTO(115)
                        !
                        ! compute the source terms
                        !
!TIMG                        call SWTSTA(117)
                        call SOURCE ( iter                , IXCGRD(1)           , IYCGRD(1)           , swpnr               , &
                                      kwave               , spcsig              , spcdir(1,2)         , spcdir(1,3)         , &
                                      ac2                 , compda(1,JDP2)      , amat(1,1,1)         , rhs                 , &
                                      abrbot              , kmespc              , dummy               , compda(1,JUBOT)     , &
                                      ufric               , compda(1,JVX2)      , compda(1,JVY2)      , idcmin              , &
                                      idcmax              , iddlow              , iddtop              , idwmin              , &
                                      idwmax              , isstop              ,                                             &
                                      swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR), &
                                      swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , swtsda(1,1,1,JP4D)  ,                       &
                                      swtsda(1,1,1,JPVEGT),                                                                   &
                                      swtsda(1,1,1,JPTRI) ,                                                                   &
                                      hs                  , etot                , qbloc               , thetaw              , &
                                      hm                  , fpm                 , wind10              , dummy               , &
                                      groww               , alimw               , smebrk              , snlc1               , &
                                      dal1                , dal2                , dal3                , ue                  , &
                                      sa1                 , sa2                 , da1c                , da1p                , &
                                      da1m                , da2c                , da2p                , da2m                , &
                                      sfnl                , dsnl                , memnl4              , wwint               , &
                                      wwawg               , wwswg               , cgo                 , compda(1,JUSTAR)    , &
                                      compda(1,JZEL)      , spcdir              , anywnd              ,                       &
                                      disc0               ,                                                                   &
                                      disc1               , genc0               , genc1               , redc0               , &
                                      redc1               , xis                 , compda(1,JFRC2)     , it                  , &
                                      compda(1,JNPLA2)    ,                                                                   &
                                      compda(1,JURSEL)    , anybin              , reflso              )
!TIMG                        call SWTSTO(117)
                        !
                     endif
                     !
                     ! update action density by means of solving the action balance equation
                     !
                     if ( ACUPDA ) then
                        !
                        ! preparatory steps before solving system of equations
                        !
!TIMG                        call SWTSTA(119)
                        call SOLPRE(ac2        , ac2old     , rhs        , amat(1,1,4), &
                                    amat(1,1,1), amat(1,1,5), amat(1,1,2), amat(1,1,3), &
                                    idcmin     , idcmax     , anybin     , idtot      , &
                                    istot      , iddlow     , iddtop     , isstop     , &
                                    spcsig     )
!TIMG                        call SWTSTO(119)
                        !
                        if ( .not.DYNDEP .and. ICUR == 0 ) then
                           !
                           ! propagation in theta space only
                           ! solve tridiagonal system of equations using Thomas' algorithm
                           !
!TIMG                           call SWTSTA(120)
                           call SOLMAT ( idcmin     , idcmax     , ac2        , rhs, &
                                         amat(1,1,1), amat(1,1,5), amat(1,1,4)     )
!TIMG                           call SWTSTO(120)
                           !
                        else
                           !
                           ! propagation in both sigma and theta spaces
                           !
                           if ( int(PNUMS(8)) == 1 ) then
                              !
                              ! implicit scheme in sigma space
                              ! solve pentadiagonal system of equations using SIP solver
                              !
!TIMG                              call SWTSTA(120)
                              call SWSIP ( ac2        , amat(1,1,1)    , rhs            , amat(1,1,4), &
                                           amat(1,1,5), amat(1,1,2)    , amat(1,1,3)    , ac2old     , &
                                           PNUMS(12)  , nint(PNUMS(14)), nint(PNUMS(13)), inocnt     , &
                                           iddlow     , iddtop         , isstop         , idcmin     , &
                                           idcmax     )
!TIMG                              call SWTSTO(120)
                              !
                           elseif (int(PNUMS(8)) == 2 ) then
                              !
                              ! explicit scheme in sigma space
                              ! solve tridiagonal system of equations using Thomas' algorithm
                              !
!TIMG                              call SWTSTA(120)
                              call SOLMT1  ( idcmin     , idcmax     , ac2        , rhs    , &
                                             amat(1,1,1), amat(1,1,5), amat(1,1,4),          &
                                             isstop     , anyblk     , iddlow     , iddtop )
!TIMG                              call SWTSTO(120)
                              !
                           endif
                           !
                        endif
                        !
                        ! if negative action density occur rescale with a factor
                        !
!TIMG                        call SWTSTA(121)
                        if ( BRESCL ) call RESCALE ( ac2, isstop, idcmin, idcmax, nrscal )
!TIMG                        call SWTSTO(121)
                        !
                        ! store propagation, generation, dissipation, redistribution, leak and radiation stress in present vertex
                        !
!TIMG                        call SWTSTA(124)
                        if ( LADDS )                                                        &
                           call ADDDIS ( compda(1,JDISS), compda(1,JLEAK), ac2   , anybin , &
                                         disc0          , disc1          ,                  &
                                         genc0          , genc1          ,                  &
                                         redc0          , redc1          ,                  &
                                         trac0          , trac1          ,                  &
                                         amat(1,1,4)    , amat(1,1,5)    ,                  &
                                         amat(1,1,2)    , amat(1,1,3)    ,                  &
                                         compda(1,JDSXB), compda(1,JDSXS), compda(1,JDSXW), &
                                         compda(1,JDSXV),                                   &
                                         compda(1,JGSXW), compda(1,JGENR),                  &
                                         compda(1,JRSXQ), compda(1,JRSXT),                  &
                                         compda(1,JREDS),                                   &
                                         compda(1,JTSXG), compda(1,JTSXT),                  &
                                         compda(1,JTSXS), compda(1,JTRAN),                  &
                                         leakcf         , compda(1,JRADS), spcsig           )
!TIMG                        call SWTSTO(124)
                        !
                        ! limit the change of the spectrum
                        !
!TIMG                        call SWTSTA(122)
                        if ( PNUMS(20) < 100. ) call PHILIM ( ac2, ac2old, cgo, kwave, spcsig, anybin, islmin, nflim, qbloc )
!TIMG                        call SWTSTO(122)
                        !
                        ! reduce the computed energy density if the value is larger then the limit value
                        ! as computed in SOURCE in case of first or second generation mode
                        !
!TIMG                        call SWTSTA(123)
                        if ( IWIND == 1 .or. IWIND == 2 ) call WINDP3 ( isstop, alimw, ac2, groww, idcmin, idcmax )
!TIMG                        call SWTSTO(123)
                        !
                        ! store some infinity norms meant for convergence check
                        !
                        if ( PNUMS(21) == 2. ) call SWACC ( ac2, ac2old, acnrms, isstop, idcmin, idcmax )
                        !
                     endif
                     !
                  endif
                  !
                  ! tag updated vertex in present cell
                  !
                  vert(ivert)%updated(jc) = icell
                  nupdv = nupdv + 1
                  !
               enddo celloop
               !
               if ( nupdv == vert(ivert)%noc ) vert(ivert)%fullupdated = .true.
!ADC               !
!ADC               ! make boundaries reflective
!ADC               if ( .false. ) then
!ADC                  !
!ADC                  if ( vert(ivert)%atti(VMARKER)==1 .and. vert(ivert)%atti(VBC)==0 ) call MakeBoundariesReflective( ivert, ac2 )
!ADC                  !
!ADC               endif
               !
            endif
            !
          enddo vertloop
          !
          ! synchronize threads before starting next sweep
          !$omp barrier
          !
       enddo sweeploop
       !
       ! global sum to inocnv counter
       !
       !$omp atomic
       inocnv = inocnv + inocnt
       !
       ! synchronize threads before checking stop criterion
       !$omp barrier
       !
!PUN       ! exchange action densities with neighbours in parallel run
!PUN       !
!PUN       if ( ACUPDA .and. MNPROC>1 ) then
!TIMG!PUN          call SWTSTA(213)
!PUN          do id = 1, MDC
!PUN             do is = 1, MSC
!PUN                ac2loc(:) = real(ac2(id,is,:),SZ)
!PUN                call UPDATER(ac2loc,dum1,dum2,1)
!PUN                ac2(id,is,:) = real(ac2loc(:),4)
!PUN             enddo
!PUN          enddo
!TIMG!PUN          call SWTSTO(213)
!PUN       endif
!PUN       !
       ! store the source terms assembled in test points per iteration in the files IFPAR, IFS1D and IFS2D
       !
       !$omp master
!TIMG       call SWTSTA(105)
       if ( NPTST > 0 .and. NSTATM == 0 ) then
          if ( IFPAR > 0 ) write (IFPAR,151) iter
          if ( IFS1D > 0 ) write (IFS1D,151) iter
          if ( IFS2D > 0 ) write (IFS2D,151) iter
          call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR), &
                        swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , &
                        swtsda(1,1,1,JPVEGT),                                                                   &
                        ac2                 , spcsig              , compda(1,JDP2)      , xytst               , &
                        dummy               )
       endif
!TIMG       call SWTSTO(105)
       !
       ! indicate number of vertices in which rescaling has taken place
       ! (only for debug purposes)
       !
       if ( ITEST >= 30 .or. idebug == 1 ) then
          !
          nwetp = 0.
          do ivert = 1, nverts
             if ( vert(ivert)%active ) nwetp = nwetp +1.
          enddo
!PUN          call SwanSumOverNodes ( nwetp )
          !
          mxnfr = maxval(nrscal)
          npfr  = count(mask=nrscal>0)
          !
!PUN          rval1 = real(npfr )
!PUN          rval2 = real(mxnfr)
!PUN          !
!PUN          ! perform global reductions in parallel run
!PUN          !
!PUN          call SwanSumOverNodes ( rval1 )
!PUN          call SwanMaxOverNodes ( rval2 )
!PUN          !
!PUN          npfr  = nint(rval1)
!PUN          mxnfr = nint(rval2)
!PUN          !
          frac = real(npfr)*100./nwetp
          if ( npfr > 0 ) write (PRINTF,134) 'rescaling', frac, mxnfr
          if ( NSTATC == 0 .and. npfr > 0 .and. IAMMASTER ) write (SCREEN,134) 'rescaling', frac, mxnfr
          !
       endif
       !
       ! indicate number of vertices in which limiter has taken place
       ! (only for debug purposes)
       !
       if ( ITEST >= 30 .or. idebug == 1 ) then
          !
          nwetp = 0.
          do ivert = 1, nverts
             if ( vert(ivert)%active ) nwetp = nwetp +1.
          enddo
!PUN          call SwanSumOverNodes ( nwetp )
          !
          mxnfl = maxval(nflim)
          npfl  = count(mask=nflim>0)
          !
!PUN          rval1 = real(npfl )
!PUN          rval2 = real(mxnfl)
!PUN          !
!PUN          ! perform global reductions in parallel run
!PUN          !
!PUN          call SwanSumOverNodes ( rval1 )
!PUN          call SwanMaxOverNodes ( rval2 )
!PUN          !
!PUN          npfl  = nint(rval1)
!PUN          mxnfl = nint(rval2)
!PUN          !
          frac = real(npfl)*100./nwetp
          if ( npfl > 0 ) write (PRINTF,134) 'limiter', frac, mxnfl
          if ( NSTATC == 0 .and. npfl > 0 .and. IAMMASTER ) write (SCREEN,134) 'limiter', frac, mxnfl
          !
          mnisl = minval(islmin)
          !
!PUN          rval1 = real(mnisl)
!PUN          call SwanMinOverNodes ( rval1 )
!PUN          mnisl = nint(rval1)
!PUN          !
          if ( npfl > 0 ) write (PRINTF,135) spcsig(mnisl)/PI2
          if ( NSTATC == 0 .and. npfl > 0 .and. IAMMASTER ) write (SCREEN,135) spcsig(mnisl)/PI2
          !
       endif
       !
       ! indicate number of vertices in which the SIP solver did not converged
       ! (only for debug purposes)
       !
       if ( ITEST >= 30 .or. idebug == 1 ) then
!PUN          rval1 = real(inocnv)
!PUN          call SwanSumOverNodes ( rval1 )
!PUN          inocnv = nint(rval1)
          if ( (DYNDEP .or. ICUR /= 0) .and. inocnv /= 0 ) then
             write (PRINTF,136) inocnv
          endif
       endif
       !$omp end master
       !
       ! info regarding the iteration process and the accuracy
       !
       if ( PNUMS(21) <= 1. ) then
          !
!TIMG          call SWTSTA(102)
          if ( PNUMS(21) == 0. ) then
             !
             call SwanConvAccur ( accur, hscurr, tmcurr, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2, ivlow, ivup )
             !
          else
             !
             call SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2, ivlow, ivup )
             !
          endif
          !
          !$omp master
          if ( iter == 1 ) then
             accur = -9999.
             write (PRINTF,142)
             if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,142)
          else
             write (PRINTF,140) accur, PNUMS(4)
             if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,140) accur, PNUMS(4)
          endif
          !$omp end master
!TIMG          call SWTSTO(102)
          !
          ! synchronize threads
          !$omp barrier
          !
          ! if accuracy has been reached then terminates iteration process
          !
          if ( accur >= PNUMS(4) ) exit iterloop
          !
       elseif ( PNUMS(21) == 2. ) then
          !
         !$omp master
!TIMG          call SWTSTA(102)
          write (PRINTF,141)
          if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,141)
          if ( iter == 1 ) then
             stopcr = -9999.
             write (PRINTF,142)
             if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,142)
          else
             if ( acnrms(1) < 1.e-20 .or. acnrmo < 1.e-20 ) then
                write (PRINTF,143)
                if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,143)
             else
                rhof   = acnrms(1)/acnrmo
                stopcr = PNUMS(1)*acnrms(2)*(1.-rhof)/rhof
                write (PRINTF,144) rhof, acnrms(1), stopcr
                if ( NSTATC == 0 .and. IAMMASTER ) write (SCREEN,144) rhof, acnrms(1), stopcr
             endif
          endif
          acnrmo = acnrms(1)
!TIMG          call SWTSTO(102)
          !$omp end master
          !
          ! synchronize threads
          !$omp barrier
          !
          ! if accuracy has been reached then terminates iteration process
          !
          if ( acnrms(1) < stopcr ) exit iterloop
          !
       endif
       !
    enddo iterloop
!TIMG    call SWTSTO(103)
    !
    ! deallocation of private arrays
    !
!TIMG    call SWTSTA(101)
    deallocate(  cad)
    deallocate(  cas)
    deallocate(  cax)
    deallocate(  cay)
    deallocate(  cgo)
    deallocate(kwave)
    !
    deallocate(idcmax)
    deallocate(idcmin)
    deallocate(iscmax)
    deallocate(iscmin)
    deallocate(anybin)
    !
    deallocate(  amat)
    deallocate(   rhs)
    deallocate(ac2old)
    !
    deallocate(anywnd)
    deallocate(obredf)
    deallocate(reflso)
    deallocate( alimw)
    deallocate( groww)
    deallocate(anyblk)
    !
    deallocate( disc0)
    deallocate( disc1)
    deallocate( genc0)
    deallocate( genc1)
    deallocate( redc0)
    deallocate( redc1)
    deallocate( trac0)
    deallocate( trac1)
    deallocate(leakcf)
    !
    deallocate(  ue)
    deallocate( sa1)
    deallocate( sa2)
    deallocate(sfnl)
    deallocate(da1c)
    deallocate(da1p)
    deallocate(da1m)
    deallocate(da2c)
    deallocate(da2p)
    deallocate(da2m)
    deallocate(dsnl)
!TIMG    call SWTSTO(101)
    !
    ! end of parallel region
    !
    !$omp end parallel
    !
    ! store the source terms assembled in test points per time step in the files IFPAR, IFS1D and IFS2D
    !
!TIMG    call SWTSTA(105)
    if ( NPTST > 0 .and. NSTATM == 1 ) then
       if ( IFPAR > 0 ) write (IFPAR,152) CHTIME
       if ( IFS1D > 0 ) write (IFS1D,152) CHTIME
       if ( IFS2D > 0 ) write (IFS2D,152) CHTIME
       call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR), &
                     swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , &
                     swtsda(1,1,1,JPVEGT),                                                                   &
                     ac2                 , spcsig              , compda(1,JDP2)      , xytst               , &
                     dummy               )
    endif
!TIMG    call SWTSTO(105)
    !
    ! deallocation of shared arrays
    !
!TIMG    call SWTSTA(101)
    deallocate(islmin)
    deallocate( nflim)
    deallocate(nrscal)
    !
    deallocate(hscurr)
    deallocate(hsprev)
    deallocate(hsdifc)
    deallocate(tmcurr)
    deallocate(tmprev)
    deallocate(tmdifc)
!PUN    deallocate(ac2loc)
    !
    deallocate(swtsda)
    !
    deallocate(memnl4)
    !
!$  deallocate(tlist)
!TIMG    call SWTSTO(101)
    !
 101 format (// ' Settings of 2nd generation mode as first guess are used:')
 102 format (// ' User-defined settings of 3rd generation mode is re-used:')
 103 format (' ITER  ',i4,'    GRWMX ',e12.4,'    ALFA   ',e12.4)
 104 format (' IWIND ',i4,'    IWCAP ',i4   ,'            IQUAD  ',i4)
 105 format (' ISURF ',i4,'    IBOT  ',i4   ,'            ITRIAD ',i4)
 106 format (' IVEG  ',i4,/)
 107 format (' Test points :',2i5)
 108 format (' Number of active points = ',i6,' (fillings-degree: ',f6.2,' %)')
 110 format ('+time ', a18, ', step ',i6, '; iteration ',i4, '; sweep ',i3)
 120 format (' iteration ', i4, '; sweep ', i3)
 130 format ('+iteration ', i4, '; sweep ', i3)
 134 format (1x,'use of ',a9,' in ',f6.2,' % of active vertices with maximum in spectral space = ',i4)
 135 format (1x,'lowest frequency occured above which limiter is applied = ',f7.4,' Hz')
 136 format (2x,'SIP solver: no convergence in ',i4,' vertices')
 140 format (' accuracy OK in ',f6.2,' % of wet grid points (',f6.2,' % required)',/ )
 141 format (7x,'Rho            ','Maxnorm   ','  Stoppingcrit.')
 142 format ('       not possible to compute accuracy, first iteration',/)
 143 format ('       norm less then 1e-20, no stopping criterion',/)
 144 format (1x,ss,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,/)
 151 format (i4, t41, 'iteration')
 152 format (a , t41, 'date-time')
    !
end subroutine SwanCompUnstruc
