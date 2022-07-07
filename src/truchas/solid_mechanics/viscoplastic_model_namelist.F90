!!
!! VISCOPLASTIC_MODEL_NAMELIST
!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_model_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_viscoplastic_model_namelists

contains

  subroutine read_viscoplastic_model_namelists(lun, params)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use truchas_logging_services
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: i_to_c
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_R

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: params

    integer :: n, ios
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: plist => null()

    !! Namelist variables
    character(128) :: phase, model
    real(r8) :: pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r
    real(r8) :: mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
        mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i
    namelist /viscoplastic_model/ phase, model, &
        pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r, &
        mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
        mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i

    call TLS_info('')
    call TLS_info('Reading VISCOPLASTIC_MODEL namelists ...')

    if (is_IOP) rewind(lun)

    n = 0
    do ! until all VISCOPLASTIC_MODEL namelists have been read or an error occurs
      if (is_IOP) call seek_to_namelist(lun, 'viscoplastic_model', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'VISCOPLASTIC_MODEL[' // i_to_c(n) // ']'

      phase = NULL_C
      model = NULL_C
      pwr_law_a = NULL_R
      pwr_law_a = NULL_R
      pwr_law_n = NULL_R
      pwr_law_q = NULL_R
      pwr_law_r = NULL_R
      mts_k = NULL_R
      mts_mu_0 = NULL_R
      mts_sig_a = NULL_R
      mts_d = NULL_R
      mts_temp_0 = NULL_R
      mts_b = NULL_R
      mts_edot_0i = NULL_R
      mts_g_0i = NULL_R
      mts_q_i = NULL_R
      mts_p_i = NULL_R
      mts_sig_i = NULL_R

      if (is_IOP) read(lun,nml=viscoplastic_model,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(phase)
      call broadcast(model)
      call broadcast(pwr_law_a)
      call broadcast(pwr_law_n)
      call broadcast(pwr_law_q)
      call broadcast(pwr_law_r)
      call broadcast(mts_k)
      call broadcast(mts_mu_0)
      call broadcast(mts_sig_a)
      call broadcast(mts_d)
      call broadcast(mts_temp_0)
      call broadcast(mts_b)
      call broadcast(mts_edot_0i)
      call broadcast(mts_g_0i)
      call broadcast(mts_q_i)
      call broadcast(mts_p_i)
      call broadcast(mts_sig_i)

      if (phase == NULL_C) then
        call TLS_fatal(label // ': PHASE not provided')
      else if (params%is_sublist(trim(phase))) then
        call TLS_fatal(label // ': another VISCOPLASTIC_MODEL namelist is already associated with this PHASE: ' // trim(phase))
      end if

      plist => params%sublist(trim(phase))
      call plist%set('model', trim(model))

      select case (trim(model))
      case ("power law")
        if (any([pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r] == NULL_R)) &
            call TLS_fatal(label // ': Missing parameters for model "power law"')
        if (any([mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
            mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i] /= NULL_R)) &
            call TLS_fatal(label // ': Given MTS parameters for model "power law"')
        call plist%set('A', pwr_law_a)
        call plist%set('n', pwr_law_n)
        call plist%set('Q', pwr_law_q)
        call plist%set('R', pwr_law_r)

      case ("MTS")
        if (any([pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r] /= NULL_R)) &
            call TLS_fatal(label // ': Given power law parameters for model "MTS"')
        if (any([mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
            mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i] == NULL_R)) &
            call TLS_fatal(label // ': Missing parameters for model "MTS"')
        call plist%set('mu0', mts_mu_0)
        call plist%set('D', mts_d)
        call plist%set('b', mts_b)
        call plist%set('k', mts_k)
        call plist%set('T0', mts_temp_0)
        call plist%set('sigma_i', mts_sig_i)
        call plist%set('sigma_a', mts_sig_a)
        call plist%set('epsdot0i', mts_edot_0i)
        call plist%set('g0i', mts_g_0i)
        call plist%set('pi', mts_p_i)
        call plist%set('qi', mts_q_i)

      case ("elastic")
        ! This is a convenience option for "ignore all parameters".
      case default
        call TLS_fatal(label // ': MODEL not recognized. Must provide "power law", "MTS", or "elastic"')
      end select
    end do

  end subroutine read_viscoplastic_model_namelists

end module viscoplastic_model_namelist
