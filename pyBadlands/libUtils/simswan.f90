!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the pyReef carbonate platform modelling application.     !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Main entry to SWAN wave modelling code.

module model

  use mpidata
  use SimWaves
  use classdata
  use miscdata

  implicit none

contains

  subroutine init(pyComm,pysFile,pysInfo,pysBot,pysOut,pyZ,pyWh,pyWp, &
                 pyWd,pyWu,pyWdd,pySide,pyDx,pyWbase,pySL,pyNx,pyNy)

      integer,intent(in) :: pyNx
      integer,intent(in) :: pyNy
      integer,intent(in) :: pyComm
      integer,intent(in) :: pySide
      real,intent(in) :: pyWu
      real,intent(in) :: pyWd
      real,intent(in) :: pyDx
      real,intent(in) :: pyWh
      real,intent(in) :: pyWp
      real,intent(in) :: pyWdd
      real,intent(in) :: pyWbase
      real,intent(in) :: pySL
      real,dimension(pyNx,pyNy),intent(in) :: pyZ
      character(len=*),intent(in) :: pysFile
      character(len=*),intent(in) :: pysInfo
      character(len=*),intent(in) :: pysBot
      character(len=*),intent(in) :: pysOut

      integer :: comm,i,j,pathlen

      ! start up MPI
      comm=pyComm
      call mpi_comm_size(comm,nprocs,ierr)
      call mpi_comm_rank(comm,iam,ierr)
      ocean_comm_world=comm

      pi=4.*atan(1.)

      ! Initialise dataset
      sp_n=pyNx
      sp_m=pyNy
      stratal_dx=pyDx
      wave_base=pyWbase
      hindcast%wvel=pyWu
      hindcast%wdir=pyWd
      hindcast%wh=pyWh
      hindcast%wp=pyWp
      hindcast%wddir=pyWdd
      sea_level = pySL

      ! if(allocated(wavU)) deallocate(wavU)
      ! allocate(wavU(pyNx,pyNy))
      ! if(allocated(wavV)) deallocate(wavV)
      ! allocate(wavV(pyNx,pyNy))

      if(allocated(sp_topo)) deallocate(sp_topo)
      allocate(sp_topo(sp_m,sp_n))
      do i=1,sp_m
         do j=1,sp_n
           sp_topo(i,j)=pyZ(j,i)-pySL
         enddo
      enddo

      ! SPM grid extent
      stratal_x=sp_n
      stratal_y=sp_m

      ! Define forecasts
      forecast_param(1)=hindcast%wh
      forecast_param(2)=hindcast%wp
      forecast_param(3)=hindcast%wdir
      forecast_param(4)=hindcast%wvel
      forecast_param(5)=hindcast%wddir
      forecast_param(6) = pySide

      int_type=mpi_integer
      real_type=mpi_real
      dbl_type=mpi_double_precision
      lgc_type=mpi_logical
      max_type=mpi_max
      min_type=mpi_min
      sum_type=mpi_sum

      pathlen=len_trim(pysFile)
      swaninput(1:pathlen) = trim(pysFile(1:pathlen))
      swaninfo = ''
      pathlen=len_trim(pysInfo)
      swaninfo(1:pathlen) = trim(pysInfo(1:pathlen))
      pathlen=len_trim(pysBot)
      swanbot = ''
      swanbot(1:pathlen) = trim(pysBot(1:pathlen))
      pathlen=len_trim(pysOut)
      swanout = ''
      swanout(1:pathlen) = trim(pysOut(1:pathlen))

      ! Initialisation step for several calls
      ! Initialise swan model dataset and directory
      call build_swan_model

      return

  end subroutine init

  subroutine run(pyComm, pyZ, pyWh, pyWp, pyWd, pyWu, pyWdd, pySide, pySL, pyWavU, &
                 pyDir, pyWavH, pyNx, pyNy)

      integer,intent(in) :: pyNx
      integer,intent(in) :: pyNy
      integer,intent(in) :: pyComm
      integer,intent(in) :: pySide
      real,intent(in) :: pyWu
      real,intent(in) :: pyWd
      real,intent(in) :: pyWh
      real,intent(in) :: pyWp
      real,intent(in) :: pyWdd
      real,intent(in) :: pySL
      real,dimension(pyNx,pyNy),intent(in) :: pyZ
      real,dimension(pyNx,pyNy),intent(out) :: pyWavU
      real,dimension(pyNx,pyNy),intent(out) :: pyDir
      real,dimension(pyNx,pyNy),intent(out) :: pyWavH
      ! real,dimension(pyNx,pyNy),intent(out) :: pyWavP
      ! real,dimension(pyNx,pyNy),intent(out) :: pyWavL

      integer :: i,j,n,p

      ocean_comm_world=pyComm

      ! Update topography
      do i=1,sp_m
         do j=1,sp_n
           sp_topo(i,j)=pyZ(j,i)-pySL
         enddo
      enddo

      ! Update forecast based on next wave regime
      forecast_param(1) = pyWh
      forecast_param(2) = pyWp
      forecast_param(3) = pyWd
      forecast_param(4) = pyWu
      forecast_param(5) = pyWdd
      forecast_param(6) = pySide

      ! Run Swan model
      call run_waves

      n=1
      do i=1,sp_m
        p=i
        do j=1,sp_n
          if(glob_Uw(p)>0.0)then
            pyWavU(j,i)=real(glob_Uw(p))
            pyDir(j,i)=real(glob_Dw(p)*pi/180.)
            !pyWavU(j,i)=real(glob_Uw(p)*cos(glob_Dw(p)*pi/180.))
            !pyWavV(j,i)=real(glob_Uw(p)*sin(glob_Dw(p)*pi/180.))
            pyWavH(j,i)=glob_Hs(p)
            !pyWavP(j,i)=glob_Per(p)
            !pyWavL(j,i)=glob_Wlen(p)
          else
            pyWavU(j,i)=0.0
            pyDir(j,i)=0.0
            pyWavH(j,i)=0.0
            !pyWavP(j,i)=0.0
            !pyWavL(j,i)=0.0
          endif
          p=p+sp_m
          n=n+1
        enddo
      enddo

      return

  end subroutine run

  subroutine final(pyComm)

      integer,intent(in) :: pyComm

      ocean_comm_world=pyComm

      ! Finalisation
      call wave_final

      return

  end subroutine final

end module model
