!                 ALLOCATABLE DATA RELATED MODULES, file 2 of 3
!
!     Contents of this file
!
!     M_WCAP             information for whitecapping
!     OUTP_DATA          information for output data
!     M_SNL4             information for quadruplets
!     M_BNDSPEC          information for boundary conditions
!     M_OBSTA            information for obstacles
!     M_GENARR           contains a number of general arrays
!     M_PARALL           information for parallelisation with MPI
!     M_DIFFR            information for diffraction
!
      MODULE M_WCAP
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
!
!  1. Updates
!
!     Sep. 00: New Module
!
!  2. Purpose
!
!     Create global variables used in whitecapping and integral parameter
!     subroutines
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
!     ACTOT  : Total action density per gridpoint
!     EDRKTOT: Zeroth moment of energy / SQRT(wavenumber)
!     EKTOT  : Zeroth moment of energy * wavenumber
!     ETOT1  : First moment of the energy density
!     ETOT2  : Second moment of the energy density
!     ETOT4  : Fourth moment of the energy density
!     KM_WAM : Mean average wavenumber according to the WAM-formulation
!     KM01   : Mean average wavenumber according to first order moment
!     SIGM_10: Mean frequency according to zeroth order moment
!     SIGM01 : Mean frequency according to first order moment
!
      REAL, SAVE    :: ACTOT
      REAL, SAVE    :: EDRKTOT
      REAL, SAVE    :: EKTOT
      REAL, SAVE    :: ETOT1
      REAL, SAVE    :: ETOT2
      REAL, SAVE    :: ETOT4
      REAL, SAVE    :: KM_WAM
      REAL, SAVE    :: KM01
      REAL, SAVE    :: SIGM_10
      REAL, SAVE    :: SIGM01
!
!$OMP THREADPRIVATE(ACTOT, EDRKTOT, EKTOT, ETOT1, ETOT2, ETOT4,
!$OMP&              KM_WAM, KM01, SIGM_10, SIGM01)
!
!     SIGPOW : contains powers of relative frequencies
!              second dimension indicates power of sigma
!
      REAL, SAVE, ALLOCATABLE :: SIGPOW(:,:)
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     SSFILL :
!     SINTGRL: Calculating integral paramters
!     WCAP   : Calculating whitecapping source term
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
      END MODULE M_WCAP

      MODULE OUTP_DATA
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
!     40.13: Nico Booij
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     40.13, July 01: New Module
!     40.13, Oct. 01: Longer filenames for output requests
!     40.31, Dec. 03: derive types OPSDAT, ORQDAT added
!
!  2. Purpose
!
!     Contains data needed during generation of output
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
      USE OCPCOMM2
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

      INTEGER, PARAMETER :: MAX_OUTP_REQ = 250 ! max. number of output requests

      CHARACTER (LEN=1)  :: OUT_COMMENT = '%' ! comment sign for heading lines

      LOGICAL            :: LCOMPGRD

      ! formats for output:
      CHARACTER (LEN=40) :: FLT_BLOCK = '(6E12.4)'       ! floating point block
      CHARACTER (LEN=40) :: FLT_TABLE = '(E11.4)'        ! floating point table
      CHARACTER (LEN=40) :: FIX_SPEC  = '(200(1X,I4))'   ! spectral output

!     format for block output per process in case of collecting data
      CHARACTER (LEN=40) :: FLT_BLKP = '(6E17.9)'

      INTEGER :: FLD_TABLE = 12       ! field length for fixed-point table
      INTEGER :: DEC_BLOCK =  4       ! number of decimals for fixed-point block
      INTEGER :: DEC_SPEC  =  4       ! number of decimals for spectral output

!     longer filenames for output requests
      CHARACTER (LEN=LENFNM) :: OUTP_FILES(1:MAX_OUTP_REQ)
      ! filenames for output; index is output request sequence number

      INTEGER, SAVE :: NREOQ = 0         ! actual number of requests saved

      TYPE OPSDAT
         CHARACTER (LEN=1)     :: PSTYPE                     ! type (F, C, P, ...)
         CHARACTER (LEN=8)     :: PSNAME                     ! name of point set
         INTEGER               :: OPI(2)                     ! integer coefficients
         REAL                  :: OPR(5)                     ! real coefficients
         INTEGER               :: MIP                        ! number of points
         REAL, POINTER         :: XP(:), YP(:), XQ(:), YQ(:) ! point coordinates
         TYPE(OPSDAT), POINTER :: NEXTOPS
      END TYPE OPSDAT

      TYPE(OPSDAT), SAVE, TARGET  :: FOPS
      TYPE(OPSDAT), SAVE, POINTER :: COPS
      LOGICAL, SAVE :: LOPS = .FALSE.

      TYPE ORQDAT
         CHARACTER (LEN=4)      :: RQTYPE   ! type (BLK, TAB, SPC ...)
         CHARACTER (LEN=8)      :: PSNAME   ! name of point set
         INTEGER                :: OQI(4)   ! integer coefficients
         REAL                   :: OQR(2)   ! real coefficients
         INTEGER, POINTER       :: IVTYP(:) ! type of output variable
         REAL, POINTER          :: FAC(:)   ! multiplication factor of block output
         TYPE(ORQDAT), POINTER  :: NEXTORQ
      END TYPE ORQDAT

      TYPE(ORQDAT), SAVE, TARGET  :: FORQ
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     SWREAD : reads data (command OUTPut OPTions)
!     SWBLOK : produces block output
!     SWTABP : produces table output
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
      END MODULE OUTP_DATA

      MODULE M_SNL4
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
!     40.17: IJsbrand Haagsma
!
!  1. Updates
!
!     Feb. 01: New Module
!
!  2. Purpose
!
!     Create global variables used in quadruplet subroutines
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

!     MDIA  : Number of quadruplets in the MDIA formulation

      INTEGER, PUBLIC, SAVE           :: MDIA = 1

!     AF11  : Contains the scaling frequency for the DIA.
!     CNL4_1: Contains the values for C1 in the MDIA formulation.
!     CNL4_2: Contains the values for C2 in the MDIA formulation.
!     LAMBDA: Contains the values for lambda in the MDIA formulation.

      REAL, PUBLIC, SAVE, ALLOCATABLE :: AF11(:)
      REAL, PUBLIC, SAVE, ALLOCATABLE :: CNL4_1(:)
      REAL, PUBLIC, SAVE, ALLOCATABLE :: CNL4_2(:)
      REAL, PUBLIC, SAVE, ALLOCATABLE :: LAMBDA(:)

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
      END MODULE M_SNL4

      MODULE M_BNDSPEC
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
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     Nov. 03: New Module
!
!  2. Purpose
!
!     Contains data with respect to specification
!     of boundary conditions
!
!  3. Method
!
!     ---
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
!     BFILED  : data concerning boundary condition files
!     BGP     : array containing data w.r.t. boundary grid points
!     BSPLOC  : place in array BSPECS where to store interpolated spectra
!     BSPDIR  : spectral directions of input spectrum
!     BSPFRQ  : spectral frequencies of input spectrum
!     CUBGP   : current item in list of boundary grid points
!     DSHAPE  : indicates option for computation of directional distribution
!               in the spectrum (boundary spectra etc.)
!               =1: directional spread in degrees is given
!               =2: power of COS is given
!     FBNDFIL : first boundary condition file in list of files
!     FBGP    : first item in list of boundary grid points
!     FBS     : first item in list of boundary spectrum parameters
!     FSHAPE  : indicates option for computation of frequency distribution
!               in the spectrum (boundary spectra etc.)
!               =1: Pierson-Moskowitz
!               =2: Jonswap
!               =3: bin
!               =4: Gaussian
!     NBS     : index of BSPECS
!     NEXTBGP : pointer to next item in list of boundary grid points
!     NEXTBS  : pointer to next item in list of boundary spectrum parameters
!     NEXTBSPC: pointer to next boundary condition file in list
!     SPPARM  : integral parameters used for computation of
!               incident spectrum. Meaning:
!               1: significant wave height
!               2: wave period (peak or mean)
!               3: average wave direction
!               4: directional distribution coefficient

      TYPE BSPCDAT
         INTEGER                :: BFILED(20)
         INTEGER, POINTER       :: BSPLOC(:)
         REAL, POINTER          :: BSPDIR(:), BSPFRQ(:)
         TYPE(BSPCDAT), POINTER :: NEXTBSPC
      END TYPE BSPCDAT

      TYPE(BSPCDAT), SAVE, TARGET :: FBNDFIL

      TYPE BSDAT
         INTEGER                :: NBS
         INTEGER                :: FSHAPE, DSHAPE
         REAL                   :: SPPARM(4)
         TYPE(BSDAT), POINTER   :: NEXTBS
      END TYPE BSDAT

      TYPE(BSDAT), SAVE, TARGET :: FBS

      TYPE BGPDAT
         INTEGER                :: BGP(6)
         TYPE(BGPDAT), POINTER  :: NEXTBGP
      END TYPE BGPDAT

      TYPE(BGPDAT), SAVE, TARGET  :: FBGP
      TYPE(BGPDAT), SAVE, POINTER :: CUBGP

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
      END MODULE M_BNDSPEC

      MODULE M_OBSTA
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
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     Nov. 03: New Module
!
!  2. Purpose
!
!     Contains data with respect to obstacles
!
!  3. Method
!
!     ---
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
!     FOBSTAC : first obstacle in list of obstacles
!     NCRPTS  : number of corner points in obstacle
!     NEXTOBST: pointer to next obstacle in list
!     RFCOEF  : reflection coefficients
!     RFTYP1  : reflection type: standard (REFL)
!     RFTYP2  : reflection type: diffusive/specular (RDIFF/RSPEC)
!     RFTYP3  : reflection type: frequency-dependent (RFD)
!     TRCOEF  : transmission coefficients
!     TRTYPE  : transmission type
!     XCRP    : x-coordinate of corner point
!     YCRP    : y-coordinate of corner point

      TYPE OBSTDAT
         INTEGER                :: TRTYPE
         REAL                   :: TRCOEF(3)
         INTEGER                :: RFTYP1, RFTYP2, RFTYP3
         REAL                   :: RFCOEF(6)
         INTEGER                :: NCRPTS
         REAL, POINTER          :: XCRP(:), YCRP(:)
         TYPE(OBSTDAT), POINTER :: NEXTOBST
      END TYPE OBSTDAT

      TYPE(OBSTDAT), SAVE, TARGET  :: FOBSTAC

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
      END MODULE M_OBSTA

      MODULE M_GENARR
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
!     40.31: Marcel Zijlema
!
!  1. Updates
!
!     Oct. 03: New Module
!
!  2. Purpose
!
!     Create several allocatable arrays for SWAN computation
!
!  3. Method
!
!     The following arrays will be created:
!
!     KGRPNT, KGRBND
!     XYTST
!     AC2
!     XCGRID, YCGRID
!     SPCSIG, SPCDIR
!     DEPTH , FRIC
!     UXB   , UYB
!     WXI   , WYI
!     WLEVL , ASTDF
!     NPLAF
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
!     AC2   : Contains action density at present time step
!     ASTDF : input field of air-sea temperature difference
!     DEPTH : input field of depth
!     FRIC  : input field of friction
!     KGRBND: array containing all boundary points
!             (+ 2 extra zeros as area separator for all separated areas)
!     KGRPNT: array containing indirect addresses for grid points
!     LAYH  : layer thickness for vegetation model
!     NPLAF : input field containing number of plants per square meter
!     SPCDIR: (*,1); spectral directions (radians)
!             (*,2); cosine of spectral directions
!             (*,3); sine of spectral directions
!             (*,4); cosine^2 of spectral directions
!             (*,5); cosine*sine of spectral directions
!             (*,6); sine^2 of spectral directions
!     SPCSIG: Relative frequencies in computational domain in sigma-space
!     UXB   : input field of contravariant U-velocity
!     UYB   : input field of contravariant V-velocity
!     VEGDIL: vegetation diameter for each layer and grid point
!     VEGDRL: drag coefficient for each layer and grid point
!     VEGNSL: number of plants / m2 for each layer and grid point
!     WLEVL : input field of water level
!     WXI   : input field of wind U-velocity (contravariant)
!     WYI   : input field of wind V-velocity (contravariant)
!     XCGRID: Coordinates of computational grid in x-direction
!     XYTST : Grid point indices of test points
!     YCGRID: Coordinates of computational grid in y-direction

      INTEGER, SAVE, ALLOCATABLE :: KGRPNT(:,:), KGRBND(:)
      INTEGER, SAVE, ALLOCATABLE :: XYTST(:)
      REAL   , SAVE, ALLOCATABLE :: AC2(:,:,:)
      REAL   , SAVE, ALLOCATABLE :: XCGRID(:,:), YCGRID(:,:)
      REAL   , SAVE, ALLOCATABLE :: SPCSIG(:)  , SPCDIR(:,:)
      REAL   , SAVE, ALLOCATABLE :: DEPTH(:) , FRIC(:)
      REAL   , SAVE, ALLOCATABLE :: UXB(:)   , UYB(:)
      REAL   , SAVE, ALLOCATABLE :: WXI(:)   , WYI(:)
      REAL   , SAVE, ALLOCATABLE :: WLEVL(:) , ASTDF(:)
      REAL   , SAVE, ALLOCATABLE :: NPLAF(:)
      REAL   , SAVE, ALLOCATABLE :: LAYH(:), VEGDIL(:), VEGDRL(:),
     &                              VEGNSL(:)

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
      END MODULE M_GENARR

      MODULE M_PARALL
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
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     Dec. 03: New Module
!     Jul. 04: introduction logicals
!
!  2. Purpose
!
!     Contains data with respect to parallel process
!     based on distributed-memory apprach using MPI
!
!  3. Method
!
!     ---
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
!     IHALOX  : width of halo area in x-direction
!     IHALOY  : width of halo area in y-direction
!     MASTER  : rank of master process
!
      INTEGER MASTER
      INTEGER IHALOX, IHALOY
      PARAMETER (MASTER = 1,
     &           IHALOX = 3, IHALOY = 3)
!
!  7. Local variables
!
!     *** variables for parallel process with MPI:
!
!     IAMMASTER: indicate whether this CPU is master or not
!     INODE   : rank of present node
!     NPROC   : number of nodes
!     PARLL   : flag to denote run as parallel (.TRUE.) or not (.FALSE.)
!     SWINT   : MPI datatype for integers
!     SWMAX   : MPI collective maximum operation
!     SWMIN   : MPI collective minimum operation
!     SWREAL  : MPI datatype for reals
!     SWSUM   : MPI collective summation
!
      INTEGER INODE, NPROC
      INTEGER SWINT, SWREAL
      INTEGER SWMAX, SWMIN, SWSUM
      LOGICAL IAMMASTER, PARLL
!
!     *** information related to global domain and subdomains
!
!     IBLKAD  : administration array for subdomain interfaces
!               contents:
!               pos. 1                     number of neighbouring subdomains
!                                          =m
!               pos. 3*i-1                 number of i-th neighbour
!               pos. 3*i                   position of i-th neighbour with
!                                          respect to present subdomain
!               pos. 3*i+1                 pointer of i-th neighbour in
!                                          last part of this array
!               pos. 3*m+2                 number of overlapping unknowns
!                                          on subdomain interface
!               pos. 3*m+3 ... 3*m+2+n     position of unknown in array
!                                          to be sent to neighbour
!               pos. 3*m+3+n ... 3*m+2*n+2 position of unknown in array
!                                          to be received from neighbour
!     IWEIG   : weights to determine load per part
!     KGRBGL  : array containing all boundary points in global domain
!               (+ 2 extra zeros as area separator for all separated areas)
!     KGRPGL  : indirect addressing for grid points in global domain
!               =1: not active point
!               >1: active point
!     LENSPO  : format length for spectral output
!     LMXF    : logical indicating whether first x-point of subdomain equals
!               first x-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMXL    : logical indicating whether last x-point of subdomain equals
!               last x-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMYF    : logical indicating whether first y-point of subdomain equals
!               first y-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMYL    : logical indicating whether last y-point of subdomain equals
!               last y-point of global domain (=.TRUE.) or not (=.FALSE.)
!     MCGRDGL : number of wet grid points in global computational grid
!     MXCGL   : number of grid points in x-direction in global
!               computational grid
!     MXF     : first index w.r.t. global grid in x-direction
!     MXL     : last index w.r.t. global grid in x-direction
!     MYCGL   : number of grid points in y-direction in global
!               computational grid
!     MYF     : first index w.r.t. global grid in y-direction
!     MYL     : last index w.r.t. global grid in y-direction
!     NBGGL   : number of grid points for which boundary condition holds
!               in global domain
!     NGRBGL  : number of boundary points in global domain
!     XCLMAX  : maximum x-coordinate in subdomain
!     XCLMIN  : minimum x-coordinate in subdomain
!     YCLMAX  : maximum y-coordinate in subdomain
!     YCLMIN  : minimum y-coordinate in subdomain
!     XGRDGL  : x-coordinate of computational grid in global domain
!     YGRDGL  : y-coordinate of computational grid in global domain
!
      INTEGER MCGRDGL, MXCGL, MYCGL
      INTEGER MXF, MXL, MYF, MYL
      INTEGER NGRBGL, NBGGL, WAV_COMM_WORLD
      REAL    XCLMAX, XCLMIN, YCLMAX, YCLMIN

      INTEGER :: LENSPO = 1000

      LOGICAL LMXF, LMXL, LMYF, LMYL

      INTEGER, SAVE, ALLOCATABLE :: IBLKAD(:)
      INTEGER, SAVE, ALLOCATABLE :: IWEIG(:)
      INTEGER, SAVE, ALLOCATABLE :: KGRPGL(:,:), KGRBGL(:)
      REAL   , SAVE, ALLOCATABLE :: XGRDGL(:,:), YGRDGL(:,:)
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
      END MODULE M_PARALL

      MODULE M_DIFFR
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
!     40.21: Agnieszka Herman, Nico Booij
!
!  1. Updates
!
!     Aug. 01: New Module
!
!  2. Purpose
!
!     Contains global variables used in diffraction procedure
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
!     DIFPARAM: contains diffraction parameter (second derivative of Hs)
!     DIFPARDX: derivative of DIFPARAM in x-direction
!     DIFPARDY: derivative of DIFPARAM in y-direction

      REAL, SAVE, ALLOCATABLE :: DIFPARAM(:)
      REAL, SAVE, ALLOCATABLE :: DIFPARDX(:)
      REAL, SAVE, ALLOCATABLE :: DIFPARDY(:)

!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     DIFPAR :   calculates the above arrays DIFPARAM, DIFPARDX, DIFPARDY
!     SPROSD :   computes propagation velocity in (x,y,theta) based on
!                arrays DIFPARAM, DIFPARDX, DIFPARDY
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
      END MODULE M_DIFFR
!/impi
!/impi      MODULE MPI
!/impi      INCLUDE 'mpif.h'
!/impi      END MODULE MPI
!PUN!/impi
!PUN!/impi      MODULE MPI
!PUN!/impi      INCLUDE 'mpif.h'
!PUN!/impi      END MODULE MPI
