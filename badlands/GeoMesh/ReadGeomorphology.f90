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
!       Filename:  ReadGeomorphology.f90
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
module readgeomorpho

  use parallel
  use FoX_sax
  use topology
  use parameters
  use hydroUtil
  use FoX_common

  implicit none

  logical,save::in_detachment=.false.
  logical,save::in_transport=.false.
  logical,save::in_lineardiff=.false.
  logical,save::in_depthdiff=.false.
  logical,save::in_nlineardiff=.false.
  logical,save::in_regolith=.false.
  logical,save::in_erodibility=.false.
  logical,save::in_perosive=.false.
  logical,save::in_spln=.false.
  logical,save::in_splm=.false.
  logical,save::in_efficiency=.false.
  logical,save::in_fracbed=.false.
  logical,save::in_stln=.false.
  logical,save::in_stlm=.false.
  logical,save::in_rock=.false.
  logical,save::in_soil=.false.
  logical,save::in_regodepth=.false.
  logical,save::in_regoprod=.false.
  logical,save::in_regothick=.false.
  logical,save::in_linear_air=.false.
  logical,save::in_linear_sea=.false.
  logical,save::in_depth_sea=.false.
  logical,save::in_depth_air=.false.
  logical,save::in_diffn=.false.
  logical,save::in_diffm=.false.
  logical,save::in_slopecritical=.false.
  logical,save::in_nlinear_sea=.false.
  logical,save::in_nlinear_air=.false.
  logical,save::in_capacity=.false.
  logical,save::in_blength=.false.
  logical,save::in_slength=.false.
  logical,save::in_chanwidth=.false.
  logical,save::in_chanexp=.false.
  logical,save::in_interface=.false.
  logical,save::in_streamero=.false.

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

    ! Detachment element
    if(name=='struct_detachment')in_detachment=.true.
    if(in_detachment) call SdetachmentElement_handler(name)
    ! Transport element
    if(name=='struct_transport')in_transport=.true.
    if(in_transport) call StransportElement_handler(name)
    ! Capacity element
    if(name=='struct_capacity')in_capacity=.true.
    if(in_capacity) call ScapacityElement_handler(name)
    ! Diffusion element
    if(name=='struct_linear_diff')in_lineardiff=.true.
    if(in_lineardiff) call SlineardiffElement_handler(name)
    if(name=='struct_depthdepend_diff')in_depthdiff=.true.
    if(in_depthdiff) call SdepthdiffElement_handler(name)
    if(name=='struct_nonlinear_diff')in_nlineardiff=.true.
    if(in_nlineardiff) call SnlineardiffElement_handler(name)
    ! Regolith element
    if(name=='struct_regolith')in_regolith=.true.
    if(in_regolith) call SregolithElement_handler(name)

  end subroutine startElement_handler
  ! =====================================================================================
  subroutine endElement_handler(namespaceURI,localname,name)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name

    call EdetachmentElement_handler(name)
    call EtransportElement_handler(name)
    call EcapacityElement_handler(name)
    call EregolithElement_handler(name)
    call ElineardiffElement_handler(name)
    call EdepthdiffElement_handler(name)
    call EnlineardiffElement_handler(name)

  end subroutine endElement_handler
  ! =====================================================================================
  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_detachment) call detachment_characters_handler(chars)
    if(in_transport) call transport_characters_handler(chars)
    if(in_capacity) call capacity_characters_handler(chars)
    if(in_regolith) call regolith_characters_handler(chars)
    if(in_lineardiff) call lineardiff_characters_handler(chars)
    if(in_depthdiff) call depthdiff_characters_handler(chars)
    if(in_nlineardiff) call nlineardiff_characters_handler(chars)

  end subroutine characters_handler
  ! =====================================================================================
  subroutine SdetachmentElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='spl_m') in_splm=.true.
    if(name=='spl_n') in_spln=.true.
    if(name=='Cerodibility') in_erodibility=.true.
    if(name=='pureErosion') in_perosive=.true.

  end subroutine SdetachmentElement_handler
  ! =====================================================================================
  subroutine EdetachmentElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='spl_m') in_splm=.false.
    if(name=='spl_n') in_spln=.false.
    if(name=='Cerodibility') in_erodibility=.false.
    if(name=='pureErosion') in_perosive=.false.

  end subroutine EdetachmentElement_handler
  ! =====================================================================================
  subroutine StransportElement_handler(name)

    character(len=*), intent(in) :: name

    if(name=='stl_m') in_stlm=.true.
    if(name=='stl_n') in_stln=.true.
    if(name=='Cefficiency') in_efficiency=.true.
    if(name=='bed_frac') in_fracbed=.true.

  end subroutine StransportElement_handler
  ! =====================================================================================
  subroutine EtransportElement_handler(name)

    character(len=*), intent(in) :: name

    if(name=='stl_m') in_stlm=.false.
    if(name=='stl_n') in_stln=.false.
    if(name=='Cefficiency') in_efficiency=.false.
    if(name=='bed_frac') in_fracbed=.false.

  end subroutine EtransportElement_handler
  ! =====================================================================================
  subroutine ScapacityElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='bedrock_length') in_blength=.true.
    if(name=='alluvial_length') in_slength=.true.
    if(name=='channel_exp') in_chanexp=.true.
    if(name=='channel_width') in_chanwidth=.true.
    if(name=='stream_erosion') in_streamero=.true.
    if(name=='bedrock_sediment_interface') in_interface=.true.

  end subroutine ScapacityElement_handler
  ! =====================================================================================
  subroutine EcapacityElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='bedrock_length') in_blength=.false.
    if(name=='alluvial_length') in_slength=.false.
    if(name=='channel_exp') in_chanexp=.false.
    if(name=='channel_width') in_chanwidth=.false.
    if(name=='stream_erosion') in_streamero=.false.
    if(name=='bedrock_sediment_interface') in_interface=.false.

  end subroutine EcapacityElement_handler
  ! =====================================================================================
  subroutine SregolithElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rego_thick') in_regothick=.true.
    if(name=='rego_prod') in_regoprod=.true.
    if(name=='rego_depth') in_regodepth=.true.
    if(name=='soil_density') in_soil=.true.
    if(name=='rock_density') in_rock=.true.

  end subroutine SregolithElement_handler
  ! =====================================================================================
  subroutine EregolithElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rego_thick') in_regothick=.false.
    if(name=='rego_prod') in_regoprod=.false.
    if(name=='rego_depth') in_regodepth=.false.
    if(name=='soil_density') in_soil=.false.
    if(name=='rock_density') in_rock=.false.

  end subroutine EregolithElement_handler
  ! =====================================================================================
  subroutine SlineardiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='coef_linear_aerial') in_linear_air=.true.
    if(name=='coef_linear_marine') in_linear_sea=.true.

  end subroutine SlineardiffElement_handler
  ! =====================================================================================
  subroutine ElineardiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='coef_linear_aerial') in_linear_air=.false.
    if(name=='coef_linear_marine') in_linear_sea=.false.

  end subroutine ElineardiffElement_handler
  ! =====================================================================================
  subroutine SdepthdiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='diff_m') in_diffm=.true.
    if(name=='diff_n') in_diffn=.true.
    if(name=='coef_depthdepend_aerial') in_depth_air=.true.
    if(name=='coef_depthdepend_marine') in_depth_sea=.true.

  end subroutine SdepthdiffElement_handler
  ! =====================================================================================
  subroutine EdepthdiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='diff_m') in_diffm=.false.
    if(name=='diff_n') in_diffn=.false.
    if(name=='coef_depthdepend_aerial') in_depth_air=.false.
    if(name=='coef_depthdepend_marine') in_depth_sea=.false.

  end subroutine EdepthdiffElement_handler
  ! =====================================================================================
  subroutine SnlineardiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='coef_non_linear_aerial') in_nlinear_air=.true.
    if(name=='coef_non_linear_marine') in_nlinear_sea=.true.
    if(name=='slope_critical') in_slopecritical=.true.

  end subroutine SnlineardiffElement_handler
  ! =====================================================================================
  subroutine EnlineardiffElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='coef_non_linear_aerial') in_nlinear_air=.false.
    if(name=='coef_non_linear_marine') in_nlinear_sea=.false.
    if(name=='slope_critical') in_slopecritical=.false.

  end subroutine EnlineardiffElement_handler
  ! =====================================================================================
  subroutine detachment_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_splm)then
      call rts(chars,spl_m)
    elseif(in_spln)then
      call rts(chars,spl_n)
    elseif(in_erodibility)then
      call rts(chars,Cerodibility)
    elseif(in_perosive)then
      call rts(chars,perosive)
    endif

  end subroutine detachment_characters_handler
  ! =====================================================================================
  subroutine transport_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_stlm)then
      call rts(chars,stl_m)
    elseif(in_stln)then
      call rts(chars,stl_n)
    elseif(in_efficiency)then
      call rts(chars,Cefficiency)
    elseif(in_fracbed)then
      call rts(chars,Fracbed)
    endif

  end subroutine transport_characters_handler
  ! =====================================================================================
  subroutine capacity_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_blength)then
      call rts(chars,bed_length)
    elseif(in_slength)then
      call rts(chars,sed_length)
    elseif(in_chanwidth)then
      call rts(chars,chan_width)
    elseif(in_chanexp)then
      call rts(chars,chan_exp)
    elseif(in_streamero)then
      call rts(chars,stream_ero)
    elseif(in_interface)then
      call rts(chars,bed_sed_interface)
    endif

  end subroutine capacity_characters_handler
  ! =====================================================================================
  subroutine regolith_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_regothick)then
      regofileflag=.true.
      regofile=chars
    elseif(in_regoprod)then
      call rts(chars,regoProd)
    elseif(in_regodepth)then
      call rts(chars,regoDepth)
    elseif(in_soil)then
      call rts(chars,soil_density)
    elseif(in_rock)then
      call rts(chars,rock_density)
    endif

  end subroutine regolith_characters_handler
  ! =====================================================================================
  subroutine lineardiff_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_linear_air)then
      call rts(chars,Cdiffusion(1))
    elseif(in_linear_sea)then
      call rts(chars,Cdiffusion(2))
    endif

  end subroutine lineardiff_characters_handler
  ! =====================================================================================
  subroutine depthdiff_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_diffm)then
      call rts(chars,Cdiff_m)
    elseif(in_diffn)then
      call rts(chars,Cdiff_n)
    elseif(in_depth_air)then
      call rts(chars,Cdiffusion_d(1))
    elseif(in_depth_sea)then
      call rts(chars,Cdiffusion_d(2))
    endif

  end subroutine depthdiff_characters_handler
  ! =====================================================================================
  subroutine nlineardiff_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_slopecritical)then
      call rts(chars,slope_critical)
    elseif(in_nlinear_air)then
      call rts(chars,Cdiffusion_nl(1))
    elseif(in_nlinear_sea)then
      call rts(chars,Cdiffusion_nl(2))
    endif

  end subroutine nlineardiff_characters_handler
  ! =====================================================================================

  subroutine geomorphology_parser

    type(xml_t)::xf

    spl_m=0.5
    spl_n=1.
    Cerodibility=0.
    perosive=0
    stl_m=1.5
    stl_n=1.
    Fracbed=0.1
    Cefficiency=0.
    Cdiffusion=0.0
    Cdiffusion_d=0.
    Cdiff_m=2.
    Cdiff_n=1.
    Cdiffusion_nl=0.0
    slope_critical=1.
    regoProd=0.
    regoDepth=0.
    soil_density=2650.
    rock_density=1325.
    bed_length=10000.
    sed_length=10000.
    stream_ero=0.
    bed_sed_interface=0.
    chan_exp=0.5
    chan_width=1.

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

  end subroutine geomorphology_parser
  ! =====================================================================================
end module readgeomorpho
! =====================================================================================
