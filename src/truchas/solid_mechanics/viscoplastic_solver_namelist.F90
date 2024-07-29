!!
!! Zach Jibben <zjibben@lanl.gov>
!! May 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_solver_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_viscoplastic_solver_namelist

contains

  subroutine read_viscoplastic_solver_namelist(lun, plist)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I, NULL_C
    use string_utilities, only: i_to_c
    use truchas_logging_services
    use physics_module, only: body_force_density

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: plist

    integer :: ios
    logical :: found
    character(80) :: iom

    real(r8) :: strain_limit ! Maximum plastic strain increment
    real(r8) :: rate_limit ! threshold for switching from Heun to BDF2
    real(r8) :: nlk_tol, nlk_vector_tolerance, relaxation_parameter, stress_relaxation_parameter
    real(r8) :: abs_plastic_strain_tol, rel_plastic_strain_tol
    integer :: maximum_iterations, nlk_max_vectors, pc_freq
    character(32) :: solver
    namelist /viscoplastic_solver/ solver, &
        strain_limit, rate_limit, &
        abs_plastic_strain_tol, rel_plastic_strain_tol, &
        nlk_tol, maximum_iterations, nlk_vector_tolerance, nlk_max_vectors, &
        pc_freq

    call TLS_info('Reading VISCOPLASTIC_SOLVER namelist ...')

    !! Locate the VISCOPLASTIC_SOLVER namelist (optional)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'viscoplastic_solver', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) then
      call TLS_info('VISCOPLASTIC_SOLVER namelist not found; using defaults.')
      return
    end if

    !! Default values
    nlk_tol = NULL_R
    maximum_iterations = NULL_I
    nlk_vector_tolerance = NULL_R
    nlk_max_vectors = NULL_I
    pc_freq = NULL_I
    relaxation_parameter = NULL_R
    stress_relaxation_parameter = NULL_R

    strain_limit = NULL_R
    rate_limit = NULL_R
    abs_plastic_strain_tol = NULL_R
    rel_plastic_strain_tol = NULL_R
    solver = NULL_C

    !! Read the namelist
    if (is_IOP) read(lun,nml=viscoplastic_solver,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading VISCOPLASTIC_SOLVER namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(maximum_iterations)
    call broadcast(nlk_tol)
    call broadcast(nlk_max_vectors)
    call broadcast(nlk_vector_tolerance)
    call broadcast(pc_freq)

    call broadcast(strain_limit)
    call broadcast(rate_limit)
    call broadcast(abs_plastic_strain_tol)
    call broadcast(rel_plastic_strain_tol)
    call broadcast(solver)

    if (strain_limit /= NULL_R) call plist%set('strain-limit', strain_limit)
    if (rate_limit /= NULL_R) call plist%set('rate-limit', rate_limit)
    if (abs_plastic_strain_tol /= NULL_R) call plist%set('atol', abs_plastic_strain_tol)
    if (rel_plastic_strain_tol /= NULL_R) call plist%set('rtol', rel_plastic_strain_tol)
    if (nlk_tol /= NULL_R) call plist%set('nlk-tol', nlk_tol)
    if (solver /= NULL_C) call plist%set('solver', solver)

    if (maximum_iterations /= NULL_I) call plist%set('nlk-max-itr', maximum_iterations)
    if (nlk_tol /= NULL_R) call plist%set('nlk-tol', nlk_tol)
    if (nlk_max_vectors /= NULL_I) call plist%set('nlk-max-vec', nlk_max_vectors)
    if (nlk_vector_tolerance /= NULL_R) call plist%set('nlk-vec-tol', nlk_vector_tolerance)
    if (pc_freq /= NULL_I) call plist%set('pc-freq', pc_freq)

  end subroutine read_viscoplastic_solver_namelist

end module viscoplastic_solver_namelist
