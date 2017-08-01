!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the pyReef carbonate platform modelling application.     !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Call to SWAN wave modelling model.

module SimWaves

    !use hdf5
    use mpidata
    use classdata
    use miscdata
    use swan_coupler_export
    use swan_coupler_parallel
    use swan_coupler_functions

    implicit none

    real::pi
    real::hs,period
    real,parameter::kappa=0.4
    real,dimension(:,:),allocatable::uCur,dCur
    real,dimension(:,:),allocatable::l_Uw,l_Dw,l_Hs,l_Per,l_Wlen
    real,dimension(:),allocatable::g_Uw,g_Dw,g_Hs,g_Per,g_Wlen
    real,dimension(:),allocatable::glob_Uw,glob_Dw,glob_Hs,glob_Per,glob_Wlen

    !  Decomposition specification for SWAN Model
    type(REEF_WADecomp),dimension(:),allocatable::wdcp

    real,dimension(6)::hcast

    public

contains

    ! ============================================================================
    subroutine build_swan_model

      integer::xpart,ypart

      call swan_createinput

      call swan_initialize(ocean_comm_world,swaninput(1:80),swaninfo(1:80))

      if(allocated(wdcp)) deallocate(wdcp)
      allocate(wdcp(nprocs))

      call swan_decomposition_grid(nprocs,wdcp)

      ! Define parameters
      xpart=wdcp(iam+1)%X_nb
      ypart=wdcp(iam+1)%Y_nb
      if(allocated(l_Uw)) deallocate(l_Uw)
      if(allocated(l_Dw)) deallocate(l_Dw)
      allocate(l_Uw(sp_n,sp_m),l_Dw(sp_n,sp_m))
      if(allocated(l_Hs)) deallocate(l_Hs)
      if(allocated(l_Wlen)) deallocate(l_Wlen)
      if(allocated(l_Per)) deallocate(l_Per)
      allocate(l_Hs(sp_n,sp_m),l_Wlen(sp_n,sp_m),l_Per(sp_n,sp_m))
      if(allocated(g_Uw)) deallocate(g_Uw)
      if(allocated(g_Dw)) deallocate(g_Dw)
      allocate(g_Uw(sp_n*sp_m),g_Dw(sp_n*sp_m))
      if(allocated(g_Hs)) deallocate(g_Hs)
      if(allocated(g_Wlen)) deallocate(g_Wlen)
      if(allocated(g_Per)) deallocate(g_Per)
      allocate(g_Hs(sp_n*sp_m),g_Wlen(sp_n*sp_m),g_Per(sp_n*sp_m))
      if(allocated(glob_Dw)) deallocate(glob_Dw)
      if(allocated(glob_Uw)) deallocate(glob_Uw)
      allocate(glob_Uw(sp_n*sp_m),glob_Dw(sp_n*sp_m))
      if(allocated(glob_Hs)) deallocate(glob_Hs)
      if(allocated(glob_Wlen)) deallocate(glob_Wlen)
      if(allocated(glob_Per)) deallocate(glob_Per)
      allocate(glob_Hs(sp_n*sp_m),glob_Wlen(sp_n*sp_m),glob_Per(sp_n*sp_m))

      return

    end subroutine build_swan_model
    ! ============================================================================
    subroutine swan_createinput

      integer::ios,iu
      real::mx_x,mx_y
      character(len=128)::stg1,stg2

      ! Create the input for SWAN simulation
      if(iam==0)then
        iu=342
        mx_x=stratal_dx*(stratal_x-1)
        mx_y=stratal_dx*(stratal_y-1)
        open(iu,file=swaninput,status="replace",action="write",iostat=ios)
        write(iu,'(a17)') "PROJECT ' ' 'S01' "
        write(iu,102) "CGRID REG",stratal_xo,stratal_yo,0,mx_x,mx_y,stratal_x-1,stratal_y-1,"CIRCLE 36 0.05 1.0"
102     format(a9,1x,f12.3,1x,f12.3,1x,i1,1x,f12.3,1x,f12.3,1x,i4,1x,i4,1x,a19)
        write(iu,103) "INPGRID BOTTOM REG",stratal_xo,stratal_yo,0,1,1,mx_x,mx_y,'EXC -999999.000'
103     format(a18,1x,f12.3,1x,f12.3,1x,i1,1x,i1,1x,i1,1x,f12.3,1x,f12.3,1x,a15)
        stg1="READINP BOTTOM 1 '"
        call append_str(stg1,swanbot)
        stg2="' 3 0 FREE"
        call append_str(stg1,stg2)
        write(iu,*)trim(stg1)
        write(iu,'(a42)') "BOUNd SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGR"

        if(forecast_param(6)==1)then
          write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE N CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==2)then
          write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE S CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==3)then
          write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE W CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==4)then
          write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE E CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==5)then
          write(iu,'(a26,3f12.3)') 'BOUNdspec SIDE NW CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==6)then
          write(iu,'(a26,3f12.3)') 'BOUNdspec SIDE NE CCW PAR ',forecast_param(1:3)
        elseif(forecast_param(6)==7)then
          write(iu,'(a26,3f12.3)') 'BOUNdspec SIDE SW CCW PAR ',forecast_param(1:3)
        else
          write(iu,'(a26,3f12.3)') 'BOUNdspec SIDE SE CCW PAR ',forecast_param(1:3)
        endif

        !write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE S CCW PAR ',forecast_param(1:3)
        !write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE N CCW PAR ',forecast_param(1:3)
        !write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE W CCW PAR ',forecast_param(1:3)
        !write(iu,'(a25,3f12.3)') 'BOUNdspec SIDE NE CCW PAR ',forecast_param(1:3)
        !write(iu,'(a7)') 'DIFFRAC'
        !write(iu,'(a8)') 'FRICTION'
        !write(iu,'(a26)') 'BREAKING CONSTANT 1.1 0.73'
        !write(iu,'(a8)') 'WCAPPING'
        !write(iu,'(a8)') 'OFF QUAD'
        write(iu,'(a4)') 'GEN1'
        write(iu,'(a26)') 'BREAKING CONSTANT 1.1 0.73'
        write(iu,'(a8)') 'FRICTION'
        !write(iu,'(a4)') 'QUAD'
        !write(iu,'(a10)') 'GEN3 AGROW'
        !write(iu,'(a10)') 'OFF BNDCHK'
        write(iu,'(a15,1x,i1,1x,i4,1x,i1,1x,i4)')"GROUP 'gf' SUBG",0,stratal_x-1,0,stratal_y-1
        stg1="TABLE 'gf' IND '"
        call append_str(stg1,swanout)
        stg2="' XP YP DIR UBOT HS PER WLEN"
        call append_str(stg1,stg2)
        write(iu,*)trim(stg1)
        write(iu,'(a5,2f12.3)') 'WIND ',forecast_param(4:5)
        write(iu,'(a7)') 'COMPUTE'
        close(iu)
        open(iu,file=swanbot,status="replace",action="write",iostat=ios)
        write(iu,'(a6)')'-1 -1'
        write(iu,'(a6)')'-1 -1'
        close(iu)
      endif
      call mpi_barrier(ocean_comm_world,ierr)

      return

    end subroutine swan_createinput
    ! ============================================================================
    subroutine run_waves

      integer::i,j,ks,ii,jj

      ! Get the bathymetry from sp model
      ii=0
      ks=iam+1
      bathyfield=0.0
      do i=wdcp(ks)%minX_id,wdcp(ks)%maxX_id
        ii=ii+1
        jj=0
        do j=wdcp(ks)%minY_id,wdcp(ks)%maxY_id
          jj=jj+1
          bathyfield(ii,jj)=-sp_topo(j,i)
          if(bathyfield(ii,jj)<0.) bathyfield(ii,jj)=-999999.0
          if(bathyfield(ii,jj)>wave_base) bathyfield(ii,jj)=-999999.0
        enddo
      enddo

      ! Define forcing waves parameters
      call import_bathymetry

      hcast(1:6)=forecast_param(1:6)
      call swan_run(hcast)
      call exportSwanData

      return

    end subroutine run_waves
    ! ============================================================================
    subroutine exportSwanData

      integer::i,j,p,n,xpart,ypart,ii,jj,xo,yo,x1,y1
      real::depth,pi

      pi=4.*atan(1.)

      xo=wdcp(iam+1)%minX_id
      yo=wdcp(iam+1)%minY_id
      x1=wdcp(iam+1)%maxX_id
      y1=wdcp(iam+1)%maxY_id
      xpart=wdcp(iam+1)%X_nb
      ypart=wdcp(iam+1)%Y_nb

      do i=1,xpart
        do j=1,ypart
          depth=bathyfield(i,j)
          ! Wave orbital velocity
          if(depth<0) exportSWANfields(i,j,1)=0.0
          if(depth<0) exportSWANfields(i,j,2)=0.0
          if(depth<0) exportSWANfields(i,j,3)=0.0
          if(depth<0) exportSWANfields(i,j,4)=0.0
          if(depth<0) exportSWANfields(i,j,5)=0.0
         enddo
      enddo
      ! Gather values to all processors
      l_Uw=-1.e8
      l_Dw=-1000.0
      l_Hs=0.
      l_Per=0.0
      l_Wlen=0.0

      ii=0
      do i=xo,x1
        ii=ii+1
        jj=0
        do j=yo,y1
          jj=jj+1
          l_Uw(i,j)=exportSWANfields(ii,jj,1)
          l_Dw(i,j)=exportSWANfields(ii,jj,2)
          l_Hs(i,j)=exportSWANfields(ii,jj,3)
          l_Per(i,j)=exportSWANfields(ii,jj,4)
          l_Wlen(i,j)=exportSWANfields(ii,jj,5)
        enddo
      enddo
      p=0
      do i=1,sp_n
        do j=1,sp_m
          p=p+1
          g_Uw(p)=l_Uw(i,j)
          g_Dw(p)=l_Dw(i,j)
          g_Hs(p)=l_Hs(i,j)
          g_Per(p)=l_Per(i,j)
          g_Wlen(p)=l_Wlen(i,j)
        enddo
      enddo
      call mpi_allreduce(g_Uw,glob_Uw,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
      call mpi_allreduce(g_Dw,glob_Dw,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
      call mpi_allreduce(g_Per,glob_Per,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
      call mpi_allreduce(g_Wlen,glob_Wlen,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
      call mpi_allreduce(g_Hs,glob_Hs,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)

      !if(iam==0) call write_h5_spm

      return

    end subroutine exportSwanData
    ! ============================================================================
    ! subroutine write_h5_spm
    !
    !   logical::compression
    !   integer::hdferr,i,j,p,n,rank
    !   character(len=128)::stg,text
    !   real(kind=8),dimension(5*sp_n*sp_m)::Vparamw
    !   integer(hid_t)::file_id
    !   integer(hid_t)::filespace,dset_id
    !   integer(hsize_t),dimension(2)::dims
    !
    !   h5data='ocean_data.h5'
    !   call noblnk(h5data)
    !   call addpath1(h5data)
    !   print*,h5data
    !   ! Initialize predefined datatypes
    !   call h5open_f(hdferr)
    !   call h5zfilter_avail_f(h5z_filter_deflate_f,compression,hdferr)
    !
    !   ! Create the file collectively.
    !   call h5fcreate_f(h5data,h5f_acc_trunc_f,file_id,hdferr)
    !
    !   dims(1)=5
    !   dims(2)=sp_n*sp_m
    !   rank=2
    !
    !   ! Create dataspace and opens it for access
    !   call h5screate_simple_f(rank,dims,filespace,hdferr)
    !   p=0
    !   do i=1,sp_m
    !     n=i
    !     do j=1,sp_n
    !       Vparamw(p+1)=glob_Uw(n)
    !       Vparamw(p+2)=glob_Dw(n)
    !       Vparamw(p+3)=glob_Hs(n)
    !       Vparamw(p+4)=glob_Per(n)
    !       Vparamw(p+5)=glob_Wlen(n)
    !       p=p+5
    !       n=n+sp_m
    !      enddo
    !   enddo
    !   text=''
    !   text="/WavesField"
    !
    !   ! Create the dataset with default properties
    !   call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,hdferr)
    !   ! Write the dataset collectively
    !   call h5dwrite_f(dset_id,h5t_native_double,Vparamw,dims,hdferr)
    !
    !   ! Close the dataset
    !   call h5dclose_f(dset_id,hdferr)
    !   call h5sclose_f(filespace,hdferr)
    !
    !   ! Close the file.
    !   call h5fclose_f(file_id,hdferr)
    !   ! Close interface
    !   call h5close_f(hdferr)
    !
    !   return
    !
    ! end subroutine write_h5_spm
    ! ============================================================================
    subroutine wave_final

      call swan_finalize
      if(allocated(bathyfield)) deallocate(bathyfield)
      if(allocated(bathyHalo)) deallocate(bathyHalo)
      if(allocated(l_Uw)) deallocate(l_Uw)
      if(allocated(g_Uw)) deallocate(g_Uw)
      if(allocated(l_Hs)) deallocate(l_Hs)
      if(allocated(g_Hs)) deallocate(g_Hs)
      if(allocated(glob_Uw)) deallocate(glob_Uw)
      if(allocated(l_Dw)) deallocate(l_Dw)
      if(allocated(g_Dw)) deallocate(g_Dw)
      if(allocated(glob_Dw)) deallocate(glob_Dw)
      if(allocated(exportSWANfields)) deallocate(exportSWANfields)
      call mpi_finalize(ierr)

      return

    end subroutine wave_final
    ! ============================================================================

end module SimWaves
