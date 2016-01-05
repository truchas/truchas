!!
!! VISCOPLASTIC_MODEL_NAMELIST
!!
!! This module provides procedures for reading the viscoplastic model
!! namelists and for accessing the data.
!!
!! Neil Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! CALL READ_VISCOPLASTIC_MODEL_NAMELISTS (LUN) reads all the so-named
!! namelists and creates the specified viscoplastic models, which are
!! held as private module data.
!!
!! GET_VP_MODEL (PHASE) returns a pointer to the viscoplastic model
!! for the specified phase.  If none was specified (meaning elastic only) a null
!! pointer is returned.  The expectation is that the solid mechanics kernel
!! uses this call only during its initialization, caching the results in its
!! own data.
!!

module viscoplastic_model_namelist

  use VP_model_class
  implicit none
  private
  
  public :: read_viscoplastic_model_namelists
  public :: get_VP_model
  
  !! The input file is read and closed as the first stage of initialization,
  !! and we need a place to hold the processed data until the solid mechanics
  !! kernel is ready to request it.  This is it; it is private data.
  !! Note that no means for deallocating this data is provided.

  type :: list_node
    integer :: seq
    character(:), allocatable :: phase
    class(VP_model), pointer :: model => null()
    type(list_node), pointer :: next  => null()
  end type
  type(list_node), pointer, save :: list => null(), last => null()
  
contains

  function get_VP_model (phase) result (model)
    character(*), intent(in) :: phase
    class(VP_model), pointer :: model
    type(list_node), pointer :: l
    l => list
    do while (associated(l))
      if (l%phase == phase) then
        model => l%model
        return
      end if
      l => l%next
    end do
    model => null()
  end function get_VP_model


  subroutine read_viscoplastic_model_namelists (lun)
  
    use kinds, only: r8
    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_C
    use parallel_communication, only: is_IOP, broadcast
    use phase_property_table
    use truchas_logging_services
    use MTS_VP_model_type
    use power_law_VP_model_type

    integer, intent(in) :: lun
    
    character(PPT_MAX_NAME_LEN) :: phase, model
    real(r8) :: pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r
    real(r8) :: mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
                mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i
    namelist /viscoplastic_model/ phase, model, pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r, &
                                  mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
                                  mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i
    
    logical :: found
    integer :: n, ios
    class(VP_model), pointer :: vpmodel
    
    call TLS_info ('')
    call TLS_info ('Reading VISCOPLASTIC_MODEL namelists ...')
    
    if (is_IOP) rewind lun
    n = 0 ! namelist counter
    
    do  ! until all VISCOPLASTICITY_MODEL namelists have been read
    
      if (is_IOP) call seek_to_namelist (lun, 'VISCOPLASTIC_MODEL', found, iostat=ios)
      call broadcast (ios)
      if (ios /= 0) call TLS_fatal ('error reading input file')
      
      call broadcast (found)
      if (.not.found) return  ! no further VISCOPLASTIC_MODEL namelists found
      
      n = n + 1
      call TLS_info ('  Reading VISCOPLASTIC_MODEL namelist #' // i_to_c(n))
      
      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        phase = NULL_C
        model = NULL_C
        pwr_law_a = NULL_R
        pwr_law_n = NULL_R
        pwr_law_q = NULL_R
        pwr_law_r = NULL_R
        mts_b = NULL_R
        mts_d = NULL_R
        mts_edot_0i = NULL_R
        mts_g_0i = NULL_R
        mts_k = NULL_R
        mts_mu_0 = NULL_R
        mts_p_i = NULL_R
        mts_q_i = NULL_R
        mts_sig_a = NULL_R
        mts_sig_i = NULL_R
        mts_temp_0 = NULL_R
        read(lun,nml=viscoplastic_model,iostat=ios)
      end if
      
      call broadcast (ios)
      if (ios /= 0) call TLS_fatal ('error reading VISCOPLASTIC_MODEL namelist')
      
      !! Broadcast the namelist variables
      call broadcast (phase)
      call broadcast (model)
      call broadcast (pwr_law_a)
      call broadcast (pwr_law_n)
      call broadcast (pwr_law_q)
      call broadcast (pwr_law_r)
      call broadcast (mts_b)
      call broadcast (mts_d)
      call broadcast (mts_edot_0i)
      call broadcast (mts_g_0i)
      call broadcast (mts_k)
      call broadcast (mts_mu_0)
      call broadcast (mts_p_i)
      call broadcast (mts_q_i)
      call broadcast (mts_sig_a)
      call broadcast (mts_sig_i)
      call broadcast (mts_temp_0)
      
      !! Check the phase name.
      if (phase == NULL_C) call TLS_fatal ('PHASE must be assigned a phase name')
      if (.not.ppt_has_phase(phase)) call TLS_fatal ('unknown PHASE: "' // trim(phase) // '"')
      if (phase_exists(phase)) &
          call TLS_fatal ('already read a VISCOPLASTIC_MODEL namelist for this PHASE: "' &
                          // trim(phase) // '"')
      
      !! Check the model name.
      if (model == NULL_C) call TLS_fatal ('MODEL must be assigned a value')
      
      select case (raise_case(model))
      case ('ELASTIC')
        vpmodel => null()
      case ('POWER LAW')
        !! Check the parameters
        if (pwr_law_a == NULL_R) call TLS_fatal ('PWR_LAW_A must be assigned a value')
        if (pwr_law_a < 0.0_r8)  call TLS_fatal ('PWR_LAW_A must be >= 0')
        if (pwr_law_n == NULL_R) call TLS_fatal ('PWR_LAW_N must be assigned a value')
        if (pwr_law_n <= 0.0_r8) call TLS_fatal ('PWR_LAW_N must be > 0')
        if (pwr_law_q == NULL_R) call TLS_fatal ('PWR_LAW_Q must be assigned a value')
        if (pwr_law_q <= 0.0_r8) call TLS_fatal ('PWR_LAW_Q must be > 0')
        if (pwr_law_r == NULL_R) call TLS_fatal ('PWR_LAW_R must be assigned a value')
        if (pwr_law_r <= 0.0_r8) call TLS_fatal ('PWR_LAW_R must be > 0')
        vpmodel => new_power_law_vp_model(pwr_law_a, pwr_law_n, pwr_law_q, pwr_law_r)
      case ('MTS')
        if (mts_b == NULL_R) call TLS_fatal ('MTS_B must be assigned a value')
        if (mts_b <= 0.0_r8) call TLS_fatal ('MTS_B must be > 0')
        if (mts_d == NULL_R) call TLS_fatal ('MTS_D must be assigned a value')
        if (mts_edot_0i == NULL_R) call TLS_fatal ('MTS_EDOT_0I must be assigned a value')
        if (mts_edot_0i < 0.0_r8) call TLS_fatal ('MTS_EDOT_0I must be >= 0')
        if (mts_g_0i == NULL_R) call TLS_fatal ('MTS_G_0I must be assigned a value')
        if (mts_g_0i <= 0.0_r8) call TLS_fatal ('MTS_G_0I must be > 0')
        if (mts_k == NULL_R) call TLS_fatal ('MTS_K must be assigned a value')
        if (mts_k <= 0.0_r8) call TLS_fatal ('MTS_K must be > 0')
        if (mts_mu_0 == NULL_R) call TLS_fatal ('MTS_MU_0 must be assigned a value')
        if (mts_mu_0 <= 0.0_r8) call TLS_fatal ('MTS_MU_0 must be > 0')
        if (mts_p_i == NULL_R) call TLS_fatal ('MTS_P_I must be assigned a value')
        if (mts_p_i <= 0.0_r8) call TLS_fatal ('MTS_P_I must be > 0')
        if (mts_q_i == NULL_R) call TLS_fatal ('MTS_Q_I must be assigned a value')
        if (mts_q_i <= 0.0_r8) call TLS_fatal ('MTS_Q_I must be > 0')
        if (mts_sig_a == NULL_R) call TLS_fatal ('MTS_SIG_A must be assigned a value')
        if (mts_sig_a < 0.0_r8) call TLS_fatal ('MTS_SIG_A must be >= 0')
        if (mts_sig_i == NULL_R) call TLS_fatal ('MTS_SIG_I must be assigned a value')
        if (mts_sig_i <= 0.0_r8) call TLS_fatal ('MTS_SIG_I must be > 0')
        if (mts_temp_0 == NULL_R) call TLS_fatal ('MTS_TEMP_0 must be assigned a value')
        if (mts_temp_0 <= 0.0_r8) call TLS_fatal ('MTS_TEMP_0 must be > 0')
        vpmodel => new_mts_vp_model(mts_k, mts_mu_0, mts_sig_a, mts_d, mts_temp_0, mts_b, &
                                    mts_edot_0i, mts_g_0i, mts_q_i, mts_p_i, mts_sig_i)
      case default
        call TLS_fatal ('Unknown viscoplastic MODEL: "' // trim(model) // '"')
      end select
      
      !! Append the data for this namelist to the list.
      if (associated(last)) then
        allocate(last%next)
        last => last%next
      else  ! list is empty
        allocate(list)
        last => list
      end if
      last%seq = n
      last%phase = trim(phase)
      last%model => vpmodel
      
    end do

  contains

    logical function phase_exists (phase)
      character(len=*), intent(in) :: phase
      type(list_node), pointer :: l
      phase_exists = .true.
      l => list
      do while (associated(l))
        if (l%phase == phase) return
        l => l%next
      end do
      phase_exists = .false.
    end function

  end subroutine read_viscoplastic_model_namelists
  
end module viscoplastic_model_namelist
