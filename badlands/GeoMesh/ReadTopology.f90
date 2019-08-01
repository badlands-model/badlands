! =====================================================================================
! BADLANDS (BAsin anD LANdscape DynamicS)
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  ReadTopology.f90
!
!    Description:  Read the XmL input file
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module readtopo

  use parallel
  use FoX_sax
  use topology
  use hydroUtil
  use parameters
  use FoX_common
  use external_forces

  implicit none

  integer::refnb,falg

  logical,save::in_outdir=.false.
  logical,save::in_geostruct=.false.
  logical,save::in_refine=.false.
  logical,save::in_regulargrid=.false.
  logical,save::in_delaunayarea=.false.
  logical,save::in_delaunayangle=.false.
  logical,save::in_refinearea=.false.
  logical,save::in_xmin=.false.
  logical,save::in_xmax=.false.
  logical,save::in_ymin=.false.
  logical,save::in_ymax=.false.
  logical,save::in_refarea=.false.
  logical,save::in_boundN=.false.
  logical,save::in_boundS=.false.
  logical,save::in_boundE=.false.
  logical,save::in_boundW=.false.
  logical,save::in_outlet=.false.
  logical,save::in_watLoss=.false.
  logical,save::in_fillAlgo=.false.
  logical,save::in_accThres=.false.
  logical,save::in_timestruct=.false.
  logical,save::in_timerestart=.false.
  logical,save::in_timefileid=.false.
  logical,save::in_procfileid=.false.
  logical,save::in_timestart=.false.
  logical,save::in_timeend=.false.
  logical,save::in_displaytime=.false.
  logical,save::in_forcestep=.false.
  logical,save::in_Tforce=.false.
  logical,save::in_udwstruct=.false.
  logical,save::in_udwtime=.false.
  logical,save::in_udwfold=.false.

  character(len=128),save::delareaS,delangleS,refineareaS
  character(len=128),save::xminS,xmaxS,yminS,ymaxS,rareaS

contains
  ! =====================================================================================
  subroutine startDocument_handler

  end subroutine startDocument_handler
  ! =====================================================================================
  subroutine endDocument_handler

  end subroutine endDocument_handler
  ! =====================================================================================
  subroutine startElement_handler(namespaceURI,localname,name,atts)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name
    type(dictionary_t),intent(in)::atts

    ! Output element
    if(name=='output_directory')in_outdir=.true.
    ! Geometry element
    if(name=='struct_geometry')in_geostruct=.true.
    if(in_geostruct) call SgeometryElement_handler(name)
    if(name=='ra') in_refine=.true.
    if(in_refine) call SrefineElement_handler(name)
    ! Time element
    if(name=='struct_time')in_timestruct=.true.
    if(in_timestruct) call StimeElement_handler(name)
    ! Underworld element
    if(name=='struct_underworld')in_udwstruct=.true.
    if(in_udwstruct) call SudwElement_handler(name)

  end subroutine startElement_handler
  ! =====================================================================================
  subroutine endElement_handler(namespaceURI,localname,name)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name

    ! Output element
    if(name=='output_directory')in_outdir=.false.
    ! Geometry element
    call EgeometryElement_handler(name)
    call ErefineElement_handler(name)
    ! Time element
    call EtimeElement_handler(name)
    ! Underworld element
    call EudwElement_handler(name)

  end subroutine endElement_handler
  ! =====================================================================================
  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    ! Get output directory name
    if(in_outdir)then
       outdir=''
       outdir=chars
    endif
    ! Geometry element
    if(in_geostruct) call geometry_characters_handler(chars)
    if(in_refine .and. refineNb>0)then
       if(.not. allocated(refine_grid)) allocate(refine_grid(refineNb))
       call refine_characters_handler(chars)
    endif
    ! Time element
    if(in_timestruct) call time_characters_handler(chars)
    ! Underworld element
    if(in_udwstruct)then
      udwFlag=.true.
      call udw_characters_handler(chars)
    endif

  end subroutine characters_handler
  ! =====================================================================================
  subroutine SgeometryElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='regular_grid') in_regulargrid=.true.
    if(name=='delaunay_area') in_delaunayarea=.true.
    if(name=='delaunay_angle') in_delaunayangle=.true.
    if(name=='refine_area') in_refinearea=.true.
    if(name=='boundN') in_boundN=.true.
    if(name=='boundS') in_boundS=.true.
    if(name=='boundE') in_boundE=.true.
    if(name=='boundW') in_boundW=.true.
    if(name=='outlet') in_outlet=.true.
    if(name=='water_loss') in_watLoss=.true.
    if(name=='fill_heigth') in_fillAlgo=.true.
    if(name=='stream_network') in_accThres=.true.

  end subroutine SgeometryElement_handler
  ! =====================================================================================
  subroutine SrefineElement_handler(name)

    character(len=*), intent(in) :: name

    if(name=='xmin') in_xmin=.true.
    if(name=='xmax') in_xmax=.true.
    if(name=='ymin') in_ymin=.true.
    if(name=='ymax') in_ymax=.true.
    if(name=='r_area') in_refarea=.true.

  end subroutine SrefineElement_handler
  ! =====================================================================================
  subroutine StimeElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='time_start') in_timestart=.true.
    if(name=='time_end') in_timeend=.true.
    if(name=='restart_folder') in_timerestart=.true.
    if(name=='restart_fileID') in_timefileid=.true.
    if(name=='restart_petNb') in_procfileid=.true.
    if(name=='display_interval') in_displaytime=.true.
    if(name=='time_step') in_forcestep=.true.
    if(name=='limit_step') in_Tforce=.true.

  end subroutine StimeElement_handler
  ! =====================================================================================
  subroutine SudwElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='sync_folder') in_udwfold=.true.
    if(name=='sync_time') in_udwtime=.true.

  end subroutine SudwElement_handler
  ! =====================================================================================
  subroutine EgeometryElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='regular_grid') in_regulargrid=.false.
    if(name=='delaunay_area') in_delaunayarea=.false.
    if(name=='delaunay_angle') in_delaunayangle=.false.
    if(name=='refine_area') in_refinearea=.false.
    if(name=='boundN') in_boundN=.false.
    if(name=='boundS') in_boundS=.false.
    if(name=='boundE') in_boundE=.false.
    if(name=='boundW') in_boundW=.false.
    if(name=='outlet') in_outlet=.false.
    if(name=='water_loss') in_watLoss=.false.
    if(name=='fill_heigth') in_fillAlgo=.false.
    if(name=='stream_network') in_accThres=.false.

  end subroutine EgeometryElement_handler
  ! =====================================================================================
  subroutine ErefineElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='xmin') in_xmin=.false.
    if(name=='xmax') in_xmax=.false.
    if(name=='ymin') in_ymin=.false.
    if(name=='ymax') in_ymax=.false.
    if(name=='r_area') in_refarea=.false.

  end subroutine ErefineElement_handler
  ! =====================================================================================
  subroutine EtimeElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='time_start') in_timestart=.false.
    if(name=='time_end') in_timeend=.false.
    if(name=='restart_folder') in_timerestart=.false.
    if(name=='restart_fileID') in_timefileid=.false.
    if(name=='restart_petNb') in_procfileid=.false.
    if(name=='display_interval') in_displaytime=.false.
    if(name=='time_step') in_forcestep=.false.
    if(name=='limit_step') in_Tforce=.false.

  end subroutine EtimeElement_handler
  ! =====================================================================================
  subroutine EudwElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='sync_folder') in_udwfold=.false.
    if(name=='sync_time') in_udwtime=.false.

  end subroutine EudwElement_handler
  ! =====================================================================================
  subroutine geometry_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_regulargrid)then
       regularfile=chars
    elseif(in_delaunayarea)then
       delareaS=chars
       call rts(delareaS,del_area)
    elseif(in_delaunayangle)then
       delangleS=chars
       call rts(delangleS,del_angle)
    elseif(in_refinearea)then
       refineareaS=chars
       call rts(refineareaS,refineNb)
    elseif(in_boundN)then
       call rts(chars,bounds(1))
    elseif(in_boundS)then
       call rts(chars,bounds(2))
    elseif(in_boundW)then
       call rts(chars,bounds(3))
    elseif(in_boundE)then
       call rts(chars,bounds(4))
    elseif(in_outlet)then
       call rts(chars,outlet)
    elseif(in_watLoss)then
       call rts(chars,infiltration_evaporation)
    elseif(in_fillAlgo)then
       call rts(chars,fh)
    elseif(in_accThres)then
      call rts(chars,accu_thres)
    endif

  end subroutine geometry_characters_handler
  ! =====================================================================================
  subroutine refine_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_xmin)then
       refnb=refnb+1
       xminS=chars
       call rts(xminS,refine_grid(refnb)%xcoord(1))
    elseif(in_xmax)then
       xmaxS=chars
       call rts(xmaxS,refine_grid(refnb)%xcoord(2))
    elseif(in_ymin)then
       yminS= chars
       call rts(yminS,refine_grid(refnb)%ycoord(1))
    elseif(in_ymax)then
       ymaxS=chars
       call rts(ymaxS,refine_grid(refnb)%ycoord(2))
    elseif(in_refarea)then
       rareaS=chars
       call rts(rareaS,refine_grid(refnb)%area)
    endif

  end subroutine refine_characters_handler
  ! =====================================================================================
  subroutine time_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_timestart)then
       call rts(chars,time_start)
    elseif(in_timeend)then
       call rts(chars,time_end)
    elseif(in_timerestart)then
        rstfolder=chars
        restartFlag=.true.
    elseif(in_timefileid)then
        call rts(chars,restartStep)
    elseif(in_procfileid)then
        call rts(chars,restartPet)
    elseif(in_displaytime)then
       call rts(chars,display_interval)
    elseif(in_forcestep)then
        call rts(chars,force_time)
    elseif(in_Tforce)then
        call rts(chars,Tforce)
    endif

  end subroutine time_characters_handler
  ! =====================================================================================
  subroutine udw_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_udwfold)then
       call rts(chars,outdir3)
    elseif(in_udwtime)then
        call rts(chars,udw_time)
    endif

  end subroutine udw_characters_handler
  ! =====================================================================================
  subroutine topology_parser

    type(xml_t)::xf
    logical::found
    integer::k
    integer::l1,l2
    character(len=128)::command,stg,fildir,file

    refnb=0

    time_start=0.
    time_end=0.
    restartPet=1
    udwFlag=.false.
    restartFlag=.false.
    display_interval=0.
    accu_thres=10
    outlet=0
    bounds=2
    fh=1.
    Tforce=0
    force_time=1.
    infiltration_evaporation=0.

    ! Default output directory
    outdir=''
    outdir(1:7)='outputs'
    outputs=''
    runfiles=''
    outputs(1:7)='outputs'
    runfiles(1:8)='runfiles'

    ! Open file
    call open_xml_file(xf,xmlfile,rc)
    if(rc/=0)then
      print*,'Failed to open XmL input file'
      call mpi_finalize(rc)
    endif

    ! Parser
    call parse(xf, &
         startDocument_handler=startDocument_handler, &
         startElement_handler=startElement_handler, &
         endElement_handler=endElement_handler, &
         characters_handler=characters_handler, &
         endDocument_handler=endDocument_handler)

    ! Close file
    call close_xml_t(xf)

    simulation_time=time_start

    ! Create simulation directory
    if(pet_id==0)then
      call noblnk(outdir)
      fildir=outdir
      stg='/runfiles/TIN.poly'
      call noblnk(stg)
      call append_str( fildir,stg )
      call noblnk(fildir)
      inquire(file=fildir,exist=found)
      if(.not.found)then
         command=' '
         command(1:6)='mkdir '
         l1=len_trim(outdir)
         command(7:l1+7)=outdir
         call term_command(command)
         command(l1+7:l1+7)='/'
         stg=''
         stg(1:l1+7)=command(1:l1+7)
         stg(l1+8:l1+15)=outputs(1:7)
         call term_command(stg)
         stg=''
         stg(1:l1+7)=command(1:l1+7)
         stg(l1+8:l1+16)=runfiles(1:8)
         call term_command(stg)
         outdir1(1:l1)=outdir(1:l1)
         outdir1(l1+1:l1+1)='/'
         outdir1(l1+2:l1+9)=outputs(1:7)
         outdir2(1:l1)=outdir(1:l1)
         outdir2(l1+1:l1+1)='/'
         outdir2(l1+2:l1+10)=runfiles(1:8)
      else
         command=' '
         command(1:6)='rm -r '
         l1=len_trim(outdir)
         command(7:l1+7)=outdir
         call term_command(command)
         command=' '
         command(1:6)='mkdir '
         l1=len_trim(outdir)
         command(7:l1+7)=outdir
         call term_command(command)
         command(l1+7:l1+7)='/'
         stg=''
         stg(1:l1+7)=command(1:l1+7)
         stg(l1+8:l1+15)=outputs(1:7)
         call term_command(stg)
         stg=''
         stg(1:l1+7)=command(1:l1+7)
         stg(l1+8:l1+16)=runfiles(1:8)
         call term_command(stg)
         outdir1(1:l1)=outdir(1:l1)
         outdir1(l1+1:l1+1)='/'
         outdir1(l1+2:l1+9)=outputs(1:7)
         outdir2(1:l1)=outdir(1:l1)
         outdir2(l1+1:l1+1)='/'
         outdir2(l1+2:l1+10)=runfiles(1:8)
      endif
      ! Underworld shared folder
      if(udwFlag)then
        call noblnk(outdir3)
        fildir=outdir3
        stg='/topsurface.vtk'
        call noblnk(stg)
        call append_str(fildir,stg)
        call noblnk(fildir)
        inquire(file=fildir,exist=found)
        if(.not.found)then
          command=' '
          command(1:6)='mkdir '
          l1=len_trim(outdir3)
          command(7:l1+7)=outdir3
          call term_command( command )
          command = ' '
          command(1:15)='chmod -R ug+rw '
          l1=len_trim(outdir3)
          command(16:l1+16)=outdir3
          call term_command( command )
        else
          command=' '
          command(1:6)='rm -r '
          l1=len_trim(outdir3)
          command(7:l1+7)=outdir3
          call term_command(command)
          command=' '
          command(1:6)='mkdir '
          l1=len_trim(outdir3)
          command(7:l1+7)=outdir3
          call term_command( command )
          command = ' '
          command(1:15)='chmod -R ug+rw '
          l1=len_trim(outdir3)
          command(16:l1+16)=outdir3
          call term_command( command )
        endif
      endif

      ! Copy XML file in the output folder
      command=''
      command(1:3)='cp '
      l1=len_trim(xmlfile)
      command(4:l1+4)=xmlfile(1:l1)
      command(l1+4:l1+4)=' '
      l2=len_trim(outdir)
      command(l1+5:l1+5+l2)=outdir(1:l2)
      command(l1+5+l2:l1+5+l2)='/'
      call term_command(command)
    endif

    call mpi_bcast(outdir1,128,mpi_character,0,badlands_world,rc)
    call mpi_bcast(outdir2,128,mpi_character,0,badlands_world,rc)

    if(udwFlag)then
      fudw='topsurface.vtk'
      fudisp='uw_output.ascii'
      maestro='maestro'
      file=''
      file=fudw
      call addpath3(file)
      fudw=file
      file=''
      file=maestro
      call addpath3(file)
      maestro=file
      file=''
      file=fudisp
      call addpath3(file)
      fudisp=file
      ! Assign displacements
      disp%event=int(time_end-time_start)
      disp%event=int(disp%event/udw_time)
      if(allocated(fgeodyn)) deallocate(fgeodyn)
      allocate(fgeodyn(disp%event))
      if(allocated(disp_time)) deallocate(disp_time)
      allocate(disp_time(disp%event,2))
      ! Broadcast vertical displacement parameter
      do k=1,disp%event
        fgeodyn(k)=fudisp
        disp_time(k,1)=time_start+udw_time*(k-1)
        disp_time(k,2)=disp_time(k,1)+udw_time
      enddo
    endif

  end subroutine topology_parser
  ! =====================================================================================
end module readtopo
! =====================================================================================
