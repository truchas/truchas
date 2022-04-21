!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_namelist

  use parameter_list_type
  implicit none
  private

  public :: read_solid_mechanics_namelist
  type(parameter_list), public :: params

contains

  subroutine read_solid_mechanics_namelist(lun)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist, NULL_R, NULL_I, NULL_C
    use string_utilities, only: i_to_c
    use truchas_logging_services
    use physics_module, only: body_force_density

    integer, intent(in) :: lun

    real(r8) :: contact_distance, contact_norm_trac, contact_penalty
    real(r8) :: strain_limit ! Maximum plastic strain increment
    real(r8) :: nlk_tol, nlk_vector_tolerance, relaxation_parameter, stress_relaxation_parameter
    real(r8) :: abs_displ_tol, rel_displ_tol, abs_stress_tol
    real(r8) :: viscoplastic_strain_limit, abs_plastic_strain_tol, rel_plastic_strain_tol, &
        viscoplastic_nlk_tol
    integer :: maximum_iterations, nlk_max_vectors, preconditioning_steps
    character(32) :: viscoplastic_solver
    namelist /solid_mechanics/ contact_distance, contact_norm_trac, contact_penalty, &
        strain_limit, nlk_tol, maximum_iterations, nlk_vector_tolerance, &
        nlk_max_vectors, preconditioning_steps, relaxation_parameter, stress_relaxation_parameter, &
        abs_displ_tol, rel_displ_tol, abs_stress_tol, &
        viscoplastic_strain_limit, abs_plastic_strain_tol, rel_plastic_strain_tol, &
        viscoplastic_nlk_tol, viscoplastic_solver

    integer :: ios
    logical :: found
    character(80) :: iom
    type(parameter_list), pointer :: plist => null()

    call TLS_info('')
    call TLS_info('Reading SOLID_MECHANICS namelist ...')

    !! Locate the SOLID_MECHANICS namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'solid_mechanics', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('SOLID_MECHANICS namelist not found')

    !! Default values
    contact_distance = NULL_R
    contact_norm_trac = NULL_R
    contact_penalty = NULL_R
    strain_limit = NULL_R

    abs_displ_tol = NULL_R
    rel_displ_tol = NULL_R
    abs_stress_tol = NULL_R

    nlk_tol = NULL_R
    maximum_iterations = NULL_I
    nlk_vector_tolerance = NULL_R
    nlk_max_vectors = NULL_I
    preconditioning_steps = NULL_I
    relaxation_parameter = NULL_R
    stress_relaxation_parameter = NULL_R

    viscoplastic_strain_limit = NULL_R
    abs_plastic_strain_tol = NULL_R
    rel_plastic_strain_tol = NULL_R
    viscoplastic_nlk_tol = NULL_R
    viscoplastic_solver = NULL_C

    !! Read the namelist
    if (is_IOP) read(lun,nml=solid_mechanics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading SOLID_MECHANICS namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call broadcast(contact_distance)
    call broadcast(contact_norm_trac)
    call broadcast(contact_penalty)
    call broadcast(strain_limit)

    call broadcast(abs_displ_tol)
    call broadcast(rel_displ_tol)
    call broadcast(abs_stress_tol)

    call broadcast(maximum_iterations)
    call broadcast(nlk_tol)
    call broadcast(nlk_max_vectors)
    call broadcast(nlk_vector_tolerance)
    call broadcast(preconditioning_steps)
    call broadcast(relaxation_parameter)
    call broadcast(stress_relaxation_parameter)

    call broadcast(viscoplastic_strain_limit)
    call broadcast(abs_plastic_strain_tol)
    call broadcast(rel_plastic_strain_tol)
    call broadcast(viscoplastic_nlk_tol)
    call broadcast(viscoplastic_solver)

    plist => params%sublist('nonlinear-solver')
    if (maximum_iterations /= NULL_I) call plist%set('nlk-max-iter', maximum_iterations)
    if (nlk_tol /= NULL_R) call plist%set('nlk-tol', nlk_tol)
    if (nlk_max_vectors /= NULL_I) call plist%set('nlk-max-vec', nlk_max_vectors)
    if (nlk_vector_tolerance /= NULL_R) call plist%set('nlk-vec-tol', nlk_vector_tolerance)

    plist => params%sublist('preconditioner')
    if (preconditioning_steps /= NULL_I) call plist%set('num-iter', preconditioning_steps)
    if (relaxation_parameter /= NULL_R) call plist%set('relaxation-parameter', relaxation_parameter)
    if (stress_relaxation_parameter /= NULL_R) &
        call plist%set('stress-relaxation-parameter', stress_relaxation_parameter)

    plist => params%sublist('model')
    if (contact_distance /= NULL_R) call plist%set('contact-distance', contact_distance)
    if (contact_norm_trac /= NULL_R) call plist%set('contact-normal-traction', contact_norm_trac)
    if (contact_penalty /= NULL_R) call plist%set('contact-penalty', contact_penalty)
    if (strain_limit /= NULL_R) call plist%set('strain-limit', strain_limit)
    if (abs_displ_tol /= NULL_R) call plist%set('abs-displ-tol', abs_displ_tol)
    if (rel_displ_tol /= NULL_R) call plist%set('rel-displ-tol', rel_displ_tol)
    if (abs_stress_tol /= NULL_R) call plist%set('abs-stress-tol', abs_stress_tol)
    call plist%set('body-force-density', body_force_density)

    plist => plist%sublist('viscoplastic-solver')
    if (viscoplastic_strain_limit /= NULL_R) call plist%set('strain-limit', viscoplastic_strain_limit)
    if (abs_plastic_strain_tol /= NULL_R) call plist%set('atol', abs_plastic_strain_tol)
    if (rel_plastic_strain_tol /= NULL_R) call plist%set('rtol', rel_plastic_strain_tol)
    if (viscoplastic_nlk_tol /= NULL_R) call plist%set('nlk-tol', viscoplastic_nlk_tol)
    if (viscoplastic_solver /= NULL_C) call plist%set('solver', viscoplastic_solver)

  end subroutine read_solid_mechanics_namelist

end module solid_mechanics_namelist
